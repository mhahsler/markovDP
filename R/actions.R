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
#' actions(Maze, state = "s(1,1)")
#' 
#' # the rest of the border can be added with more entries in the reward 
#' # function. But since the algorithm learns not to waste moves, the policy
#' # eventually will not go out of boundary anyway.
#' @returns a vector with the available actions.
#' @export
actions <- function(x, state) {
  a <- x$actions[!sapply(x$actions, FUN = function(a) {
    all(reward_matrix(x, action = a, start.state = state) == -Inf)
  })]
  
  a <- factor(a, levels = x$actions)
  a
}