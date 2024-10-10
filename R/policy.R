#' Extract or Create a Policy
#'
#' Extracts the policy from a solved model or create a policy. All
#' policies are deterministic.
#'
#' For an MDP, the deterministic policy is a data.frame with columns for:
#'
#' * `state`: The state.
#' * `U`: The state's value (discounted expected utility U) if the policy
#'    is followed.
#' * `action`: The prescribed action.
#'
#' For unconverged, finite-horizon problems, the solution is a policy for
#' each epoch. This is returned as a list of data.frames.
#'
#' @family policy
#'
#' @param model A solved [MDP] object.
#' @param epoch return the policy of the given epoch. `NULL` returns a list
#'   with elements for each epoch.
#' @param drop logical; drop the list for converged, epoch-independent policies.
#' @return A data.frame containing the policy. If `drop = FALSE` then the policy is returned
#'  as a list with the policy for each epoch.
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
#' ## create a random policy
#' pi_random <- random_policy(Maze, estimate_U = TRUE)
#' pi_random
#'
#' gw_plot(add_policy(Maze, pi_random))
#'
#' ## create a manual policy (go up and in some squares to the right)
#' acts <- rep("up", times = length(Maze$states))
#' names(acts) <- Maze$states
#' acts[c("s(1,1)", "s(1,2)", "s(1,3)")] <- "right"
#' pi_manual <- manual_policy(Maze, acts)
#' pi_manual
#'
#' gw_plot(add_policy(Maze, pi_manual))
#'
#' ## Finite horizon (we use incremental pruning because grid does not converge)
#' sol <- solve_MDP(model = Maze, horizon = 3)
#' sol
#'
#' policy(sol)
#' gw_plot(sol)
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
#' @param prob probability vector for random actions for `random_policy()`.
#'   a logical indicating if action probabilities should be returned for
#'   `greedy_action()`.
#' @param estimate_U logical; estimate the value function
#'    using [policy_evaluation()]?
#' @param only_available_actions logical; only sample from available actions?
#'   (see [available_actions()] for details)
#' @param ... is passed on to [available_actions()].
#' @export
random_policy <-
  function(model,
           prob = NULL,
           estimate_U = FALSE,
           only_available_actions = FALSE,
           ...) {
    if (!inherits(model, "MDP")) {
      stop("'model' needs to be of class 'MDP'.")
    }
    
    A <- model$actions
    S <- model$states
    
    if (only_available_actions) {
      if (is.null(prob))
        prob <- rep(1/length(A), times = length(A))
      
      pol <- data.frame(
        state = S,
        U = NA_real_,
        action = .normalize_action(
          sapply(S, function(s) {
            avail_a <- available_actions(model, s, ...)
            if (length(avail_a) == 0L)
              avail_a <- seq_along(A)
            if (length(avail_a) == 1L)
              avail_a
            else 
              sample(as.integer(avail_a), size = 1L, prob = prob[avail_a]/sum(prob[avail_a]))
          }), model
        ), row.names = NULL
      )
      
      
    } else {
      pol <- data.frame(
        state = S,
        U = NA_real_,
        action = .normalize_action(
          sample(
            seq_along(A),
            size = length(S),
            replace = TRUE,
            prob = prob
          ), model
        ), row.names = NULL
      )
    }
    
    if (estimate_U) {
      pol$U <- policy_evaluation(model, pol)
    }
    
    pol
  }


#' @rdname policy
#' @param actions a vector with the action (either the action label or the
#'  numeric id) for each state.
#' @param U a vector with the value function for the policy. If `TRUE`, then
#'    the it is estimated using `policy_evaluation()`.
#' @export
manual_policy <-
  function(model, actions, U = NULL, estimate_U = FALSE) {
    if (!inherits(model, "MDP")) {
      stop("'model' needs to be of class 'MDP'.")
    }
    
    if (is.null(U)) {
      U <- NA_real_
    } else {
      estimate_U <- FALSE
    }
    
    A <- model$actions
    S <- model$states
    
    if (is.numeric(actions)) {
      actions <- A[actions]
    }
    
    pol <- data.frame(state = S,
                      U = U,
                      action = .normalize_action(actions, model), 
                      row.names = NULL)
    
    if (estimate_U) {
      pol$U <- policy_evaluation(model, pol)
    }
    
    pol
  }
