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
#' gridworld_plot(sol)
#'
#' ## create a random policy
#' pi_random <- random_policy(Maze)
#' pi_random
#'
#' gridworld_plot(add_policy(Maze, pi_random))
#'
#' ## create a manual policy (go up and in some squares to the right)
#' acts <- rep("up", times = length(Maze$states))
#' names(acts) <- Maze$states
#' acts[c("s(1,1)", "s(1,2)", "s(1,3)")] <- "right"
#' pi_manual <- manual_policy(Maze, acts)
#' pi_manual
#'
#' gridworld_plot(add_policy(Maze, pi_manual))
#'
#' ## Finite horizon (we use incremental pruning because grid does not converge)
#' sol <- solve_MDP(model = Maze, horizon = 3)
#' sol
#'
#' policy(sol)
#' gridworld_plot(sol)
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
#' @param prob probability vector for random actions for `random_policy()`.
#'   a logical indicating if action probabilities should be returned for
#'   `greedy_action()`.
#' @export
random_policy <-
  function(model, prob = NULL, U = FALSE) {
    if (!inherits(model, "MDP")) {
      stop("'model' needs to be of class 'MDP'.")
    }
    
    A <- model$actions
    S <- model$states
    
    pol <- data.frame(
      state = S,
      U = NA_real_,
      action = factor(
        sample(
          seq_along(A),
          size = length(S),
          replace = TRUE,
          prob = prob
        ),
        levels = seq_along(A),
        labels = A
      )
    )
    
    if (U) {
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
  function(model, actions, U = NULL) {
    if (!inherits(model, "MDP")) {
      stop("'model' needs to be of class 'MDP'.")
    }
    
    estimate_U <- FALSE
    if (is.null(U))
      U <- NA_real_
    else if (is.logical(U)) {
      if(U) 
        estimate_U <- TRUE
      U <- NA_real_  
    }
    
    A <- model$actions
    S <- model$states
    
    if (is.numeric(actions)) {
      actions <- A[actions]
    }
    
    actions <- factor(actions, levels = A)
    
    
    pol <- data.frame(
      state = S,
      U = U,
      action = actions
    )
    
    if (all(is.na(U)) && estimate_U) {
      pol$U <- policy_evaluation(model, pol)
    }
    
    pol
  }
