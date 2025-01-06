#' Extract, Create Add a Policy to a Model
#'
#' Extracts the policy from a solved model or create a policy. All
#' policies are deterministic.
#'
#' `policy()` extracts the (deterministic) policy from a solved MDP in the form of a 
#' a data.frame with columns for:
#'
#' * `state`: The state.
#' * `V`: The state values if the policy is followed.
#' * `action`: The prescribed action.
#'
#' For unconverged, finite-horizon problems, the solution is a policy for
#' each epoch. This is returned as a list of data.frames.
#'
#' `add_policy()` adds a policy to an existing MDP object.
#'
#' `random_policy()` and `manual_policy()` construct new policies.
#' 
#' `induced_transition_matrix()` returns the single transition matrix which follows
#' the actions specified in a policy.
#' @family policy
#'
#' @param model A solved [MDP] object.
#' @param epoch return the policy of the given epoch. `NULL` returns a list
#'   with elements for each epoch.
#' @param drop logical; drop the list for converged, epoch-independent policies.
#' @returns 
#'  - `policy()`, `random_policy()` and `manual_policy()` return a data.frame 
#'    containing the policy. If `drop = FALSE` then the policy is returned
#'    as a list with the policy for each epoch.
#'  - `add_policy()` returns an MDP object.
#'  - `induced_transition_matrix` returns a single transition matrix.
#' @author Michael Hahsler
#' @keywords graphs
#' @examples
#' data("Maze")
#'
#' sol <- solve_MDP(Maze)
#' sol
#'
#' ## policy with value function and optimal action.
#' policy(sol)
#' plot_value_function(sol)
#' gw_plot(sol)
#'
#' induced_transition_matrix(sol)
#'
#' ## create a random policy
#' pi_random <- random_policy(Maze, estimate_V = TRUE)
#' pi_random
#'
#' gw_plot(add_policy(Maze, pi_random))
#'
#' ## create a manual policy (go up and in some squares to the right)
#' acts <- rep("up", times = length(Maze$states))
#' names(acts) <- Maze$states
#' acts[c("s(1,1)", "s(1,2)", "s(1,3)")] <- "right"
#' acts
#' 
#' pi_manual <- manual_policy(Maze, acts, estimate_V = TRUE)
#' pi_manual
#'
#' gw_plot(add_policy(Maze, pi_manual))
#'
#' # Transition matrix induced by the policy
#' induced_transition_matrix(Maze, pi_manual, sparse = TRUE)
#'
#' ## Finite horizon (we use incremental pruning because grid does not converge)
#' sol <- solve_MDP(model = Maze, horizon = 3)
#' sol
#'
#' policy(sol)
#' gw_plot(sol, epoch = 1)
#' gw_plot(sol, epoch = 2)
#' gw_plot(sol, epoch = 3)
#' @export
policy <- function(model, epoch = NULL, drop = TRUE) {
  UseMethod("policy")
}

#' @export
policy.MDP <- function(model, epoch = NULL, drop = TRUE) {
  is_solved_MDP(model, stop = TRUE)
  
  policy <- model$solution$policy
  
  if (!is.null(epoch)) {
    return(policy[[.get_pol_index(model, epoch)]])
  }
  
  if (drop && length(policy) == 1) {
    policy <- policy[[1]]
  }
  
  policy
}

#' @rdname policy
#'
#' @param policy a policy data.frame.
#'
#' @returns The model description with the added policy.
#'
#' @author Michael Hahsler
#' @export
add_policy <- function(model, policy) {
  UseMethod("add_policy")
}

#' @export
add_policy.MDP <- function(model, policy) {
  policy <- .normalize_policy(policy, model)
  
  solution <- list(method = "manual",
                   policy = list(policy),
                   converged = NA)
  
  model$solution <- solution
  #model <- check_and_fix_MDP(model)
  
  model
}

.normalize_policy <- function(policy,
                              model,
                              V = NULL,
                              data_frame = TRUE) {
  if (!is.null(V) && !((length(V) == 1L) || length(V) == length(S(model))))
    stop("Specified value function is not compartible with the number of states!")
  
  if (is.data.frame(policy)) {
    action <- .normalize_action(policy$action, model)
    V <- V %||% policy$V
  }  else {
    action <- .normalize_action(policy, model)
    V <- V %||% NA_real_
  }
  
  # checks
  if (length(action) != length(S(model)))
    stop("Policy length does not match number of states!")
  
  if (any(is.na(action)))
    stop("Actions do not match the model!")
  
  # TODO policy eval?
  if (data_frame)
    return(data.frame(
      states = S(model),
      V = V %||% NA_real_,
      action = action,
      row.names = seq_along(S(model))
    ))
  else
    return(action)
}

#' @rdname policy
#'
#' @param prob probability vector for random actions for `random_policy()`.
#'   a logical indicating if action probabilities should be returned for
#'   `greedy_action()`.
#' @param estimate_V logical; estimate the value function
#'    using [policy_evaluation()]?
#' @param only_available_actions logical; only sample from available actions?
#'   (see [available_actions()] for details)
#' @param ... is passed on to [available_actions()].
#' @export
random_policy <-
  function(model,
           prob = NULL,
           estimate_V = FALSE,
           only_available_actions = FALSE,
           ...) {
    if (!inherits(model, "MDP")) {
      stop("'model' needs to be of class 'MDP'.")
    }
    
    A <- A(model)
    S <- S(model)
    
    if (only_available_actions) {
      if (is.null(prob))
        prob <- rep(1 / length(A), times = length(A))
      
      action <- sapply(S, function(s) {
        avail_a <- available_actions(model, s, ...)
        if (length(avail_a) == 0L)
          avail_a <- seq_along(A)
        if (length(avail_a) == 1L)
          avail_a
        else
          sample(as.integer(avail_a),
                 size = 1L,
                 prob = prob[avail_a] / sum(prob[avail_a]))
      })
    } else {
      action <- sample.int(length(A),
                           size = length(S),
                           replace = TRUE,
                           prob = prob)
    }
    
    policy <- .normalize_policy(action, model, data_frame = TRUE)
    
    if (estimate_V) {
      policy$V <- policy_evaluation(model, policy)
    }
    
    policy
  }


#' @rdname policy
#' @param actions a vector with the action (either the action label or the
#'  numeric id) for each state.
#' @param V a vector representing the value function for the policy. If `TRUE`, then
#'    the it is estimated using `policy_evaluation()`.
#' @param sparse logical; should a sparse transition matrix be returned?
#' @export
manual_policy <-
  function(model,
           actions,
           V = NULL,
           estimate_V = FALSE) {
    if (!inherits(model, "MDP")) {
      stop("'model' needs to be of class 'MDP'.")
    }
    
    if (is.null(V)) {
      V <- NA_real_
    } else {
      estimate_V <- FALSE
    }
    
    if (is.numeric(actions)) {
      actions <- A(model)[actions]
    }
    
    policy <- .normalize_policy(actions, model, data_frame = TRUE)
    
    if (estimate_V) {
      policy$V <- policy_evaluation(model, policy)
    }
    
    policy
  }

#' @rdname policy
#' @export
induced_transition_matrix <- function(model, policy = NULL, epoch = 1L, sparse = FALSE) {
  if (is.null(policy))
    policy <- policy(model, epoch = epoch)
  
  if (is.data.frame(policy))
    policy <- policy$action
    
  P <- sapply(seq_along(S(model)), FUN = function(s) transition_matrix(model, policy[s], s))
  if (is.list(P) && is(P[[1]], "sparseVector"))
    P <- do.call(MatrixExtra::rbind_csr, P)
  else
    P <- t(P)
  dimnames(P) <- list(S(model), S(model))
 
  P <- .sparsify(P, sparse)
  P
}

#' @rdname policy
#' @export
induced_reward_matrix <-
  function(model, policy = NULL, epoch = 1L) {
    if (is.null(policy))
      policy <- policy(model, epoch = epoch)
    
    if (is.data.frame(policy))
      policy <- policy$action
    
    t(sapply(seq_along(S(model)), FUN = function(s) reward_matrix(model, policy[s], s)))
  }
