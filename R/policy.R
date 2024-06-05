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
#' @param x A solved [MDP] object.
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
policy <- function(x, epoch = NULL, drop = TRUE) {
  UseMethod("policy")
}

#' @export
policy.MDP <- function(x, epoch = NULL, drop = TRUE) {
  is_solved_MDP(x, stop = TRUE)

  policy <- x$solution$policy

  if (!is.null(epoch)) {
    return(policy[[.get_pol_index(x, epoch)]])
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
  function(x, prob = NULL) {
    if (!inherits(x, "MDP")) {
      stop("'x' needs to be of class 'MDP'.")
    }
    
    A <- x$actions
    S <- x$states
    
    data.frame(
      state = S,
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
  }


#' @rdname policy
#' @param actions a vector with the action (either the action label or the
#'  numeric id) for each state.
#' @export
manual_policy <-
  function(x, actions) {
    if (!inherits(x, "MDP")) {
      stop("'x' needs to be of class 'MDP'.")
    }
    
    A <- x$actions
    S <- x$states
    
    if (is.numeric(actions)) {
      actions <- A[actions]
    }
    
    actions <- factor(actions, levels = A)
    
    data.frame(
      state = S,
      action = actions
    )
  }
