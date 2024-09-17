#' Available Actions in a State
#'
#' Determine the set of actions available in a state.
#'
#' Unavailable actions are modeled as actions that have an immediate
#'  reward of `-Inf` in the reward function.
#' @name actions
#' @family MDP
#'
#' @param x a [MDP] object.
#' @param state a character vector of length one specifying the state.
#' @param inf_reward logical; consider an action that produced `-Inf` reward to all 
#'  end states unavailable?
#' @param stay_in_place logical; consider an action that results in the same state
#'  with a probability of 1 as unavailable. Note that this will mean that
#'  absorbing states have not available action!
#' @returns a character vector with the available actions.
#'
#' @author Michael Hahsler
#' @examples
#' data(Maze)
#' gridworld_matrix(Maze)
#'
#' # The the following actions are always available:
#' Maze$actions
#'
#' # An action that leaves the grid currently is allowed but does not do 
#' # anything.
#' act(Maze, "s(1,1)", "up")
#'
#' # Make the action unavailable by setting the reward to -Inf
#' Maze$reward
#' Maze$reward <- rbind(
#'     Maze$reward, 
#'     R_(action = "up", start.state = "s(1,1)", value = - Inf))
#'
#' # up in s(1,1) now produces a reward of - Inf
#' act(Maze, "s(1,1)", "up")
#'
#' # up is unavailable for s(1,1)
#' available_actions(Maze, state = "s(1,1)")
#' 
#' # the rest of the border can be added with more entries in the reward 
#' # function. But since the algorithm learns not to waste moves, the policy
#' # eventually will not go out of boundary anyway.
#' @returns a vector with the available actions.
#' @export
available_actions <- function(x, state, inf_reward = TRUE, stay_in_place = TRUE) {
 
  acts <- x$actions
  
  if (stay_in_place) {
    acts <- acts[transition_matrix(
      x,
      action = acts,
      start.state = state,
      end.state = state,
      simplify = TRUE
    ) != 1L]
  }
   
  if (inf_reward) {
    acts <- acts[!sapply(
      acts,
      FUN = function(a) {
        all(reward_matrix(x, action = a, start.state = state) == -Inf)
      }
    )]
  }
  
  acts <- factor(acts, levels = x$actions)
  acts
}