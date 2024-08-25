#' Perform an Action
#'
#' Performs an action in a state and returns the new state and reward.
#'
#' @family MDP
#'
#' @param model a MDP model.
#' @param state the current state.
#' @param action the chosen action.
#' 
#' @return a names list with the next `state` and the `reward`.
#'
#' @author Michael Hahsler
#' @examples
#' data(Maze)
#'
#' act(Maze, "s(1,3)", "right")
#' @export
act <- function(model, state, action) {
  S <- model$states
  sp <- sample.int(length(S), 
                   1L, 
                   prob = transition_matrix(model, action, state))
  sp <- factor(S[sp], levels = S)
  r <-  reward_matrix(model, action, state, sp)
                   
  list(state = sp,
       reward =  r)
}

