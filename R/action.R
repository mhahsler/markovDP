#' Action Given a Policy
#'
#' Returns the action given a policy. If the 
#' policy is optimal, then also the action will be optimal.
#'
#' @family policy

#' @param model a solved [MDP].
#' @param state the state.
#' @param epoch what epoch of the policy should be used. Use 1 for converged policies.
#' @param ... further parameters are passed on.
#' @return The name of the optimal action.
#' @author Michael Hahsler
#' @examples
#' data("Maze")
#' Maze
#'
#' sol <- solve_MDP(Maze)
#' policy(sol)
#'
#' action(sol, state = "s(1,3)")
#' @export
action <- function(model, ...) {
  UseMethod("action")
}

#' @rdname action
#' @export
action.MDP <-
  function(model, state, epoch = 1, ...) {
    policy(model, epoch)$action[.get_state_id(model, state)]
}