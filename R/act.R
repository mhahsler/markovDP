#' Perform an Action
#'
#' Performs an action in a state and returns the new state and reward.
#'
#' @family MDP
#'
#' @param model an MDP model.
#' @param state the current state.
#' @param action the chosen action. If the action is not specified (`NULL`) and
#'   the MDP model contains a policy, then the action is chosen according to the
#'   policy.
#' @param ... if action is unspecified, then the additional parameters are 
#'   passed on to `action()` to determine the action using the model's policy.
#' 
#' @returns a names list with the `old_state`, the `action`, 
#'  the next `reward` and `state`.
#'
#' @author Michael Hahsler
#' @examples
#' data(Maze)
#'
#' act(Maze, "s(1,3)", "right")
#' 
#' # solve the maze and then ask for actions using the policy
#' sol <- solve_MDP(Maze)
#' act(sol, "s(1,3)")
#' 
#' # make the policy in sol epsilon-soft and ask 10 times for the action 
#' replicate(10, act(sol, "s(1,3)", epsilon = .2))
#' @export
act <- function(model, state, action = NULL, ...) {
  # follow the policy in the model?
  action <- action %||% action(model, state, ...)
  
  action <- .normalize_action(action, model)
  state <- .normalize_state(state, model)
  
  sp <- sample.int(length(S(model)), 
                   1L, 
                   prob = transition_matrix(model, action, state, 
                                            sparse = FALSE))
  
  sp <- .normalize_state(sp, model)
  r <-  reward_matrix(model, action, state, sp)
  
  list(old_state = state,
       action = action,
       reward =  r,
       state = sp)
}

# faster without conversions
.act_int <- function(model, state, action, ...) {
  sp <- sample.int(length(S(model)), 
                   1L, 
                   prob = transition_matrix(model, action, state, 
                                            sparse = FALSE))

  r <-  reward_matrix(model, action, state, sp)
  list(old_state = state,
       action = action,
       reward =  r,
       state = sp)
}
