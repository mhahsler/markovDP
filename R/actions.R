#' Available Actions in a State
#'
#' Determine the set of actions available in a state.
#'
#' Unavailable actions are modeled here a actions that have an immediate
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
#' # available actions in
#' actions(Maze, state = "s(3,1)")
#' @returns a vector with the available actions.
#' @export
actions <- function(x, state) {
  a <- x$actions[!sapply(x$actions, FUN = function(a) {
    all(reward_matrix(x, action = a, start.state = state) == -Inf)
  })]
  
  a <- factor(a, levels = x$actions)
  a
}