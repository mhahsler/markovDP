#' Extract the Policy 
#'
#' Extracts the policy from a solved model.
#'
#' A list (one entry per epoch) with the optimal policy.
#' For converged, infinite-horizon problems solutions, a list with only the
#' converged solution is produced.
#' 
#' For an MDP, the policy is a data.frame with columns for:
#'
#' * `state`: The state.
#' * `U`: The state's value (discounted expected utility U) if the policy
#'    is followed
#' * `action`: The prescribed action.
#'
#' @family policy
#'
#' @param x A solved [MDP] object.
#' @param epoch return the policy of the given epoch. `NULL` returns a list 
#'   with elements for each epoch.
#' @param drop logical; drop the list for converged, epoch-independent policies.
#' @return A list with the policy for each epoch. Converged policies
#'  have only one element. If `drop = TRUE` then the policy is returned
#'  without a list.
#' @author Michael Hahsler
#' @keywords graphs
#' @examples
#' data("Maze")
#'
#' sol <- solve_MDP(Maze)
#' sol
#'
#' # policy with value function and optimal action.
#' policy(sol)
#' plot_value_function(sol)
#'
#' # Finite horizon (we use incremental pruning because grid does not converge)
#' sol <- solve_MDP(model = Maze, horizon = 3)
#' sol
#'
#' policy(sol)
#' @export
policy <- function(x, epoch = NULL, drop = TRUE) {
  UseMethod("policy")
}

#' @export
policy.MDP <- function(x, epoch = NULL, drop = TRUE) {
  is_solved_MDP(x, stop = TRUE)

  policy <- x$solution$policy
  
  if (!is.null(epoch))
    return (policy[[.get_pol_index(x, epoch)]])
  
  if (drop && length(policy) == 1) {
    policy <- policy[[1]]
  }

  policy
}
