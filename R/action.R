#' Action Given a Policy
#'
#' Returns an action given a policy. If the
#' policy is optimal, then also the action will be optimal.
#'
#' @family policy

#' @param model a solved [MDP].
#' @param state the state.
#' @param epoch what epoch of the policy should be used. Use 1 for converged policies.
#' @param epsilon make the policy epsilon soft.
#' @param ... further parameters are passed on.
#' @return The name of the optimal action as a factor.
#' @author Michael Hahsler
#' @examples
#' data("Maze")
#' Maze
#'
#' sol <- solve_MDP(Maze)
#' policy(sol)
#'
#' action(sol, state = "s(1,3)")
#' 
#' ## choose from an epsilon-soft policy
#' table(replicate(100, action(sol, state = "s(1,3)", epsilon = 0.1)))
#' @export
action <- function(model, ...) {
  UseMethod("action")
}

#' @rdname action
#' @export
action.MDP <-
  function(model, state, epsilon = 0, epoch = 1, ...) {
    a <- policy(model, epoch)$action[.get_state_id(model, state)]
    
    if (epsilon == 0)
      return(a)
    
    available_A <- actions(model, state) 
      
    if (length(available_A) > 1L && runif(1) < epsilon) {
      a <- sample(available_A, size = 1L)
    }
   
    a 
  }
