#' Calculate the Expected Reward of a Policy
#'
#' This function calculates the expected total reward for a MDP policy
#' given a start state (distribution). The value is calculated using the value function stored
#' in the MDP solution.
#'
#' The reward is typically calculated using the value function
#' of the solution. If these are not available, then [simulate_MDP()] is
#' used instead with a warning.
#'
#' @family policy
#'
#' @param x a solved [MDP] object.
#' @param start specification of the current state (see argument start
#' in [MDP] for details). By default the start state defined in
#' the model as start is used. Multiple states can be specified as rows in a matrix.
#' @param epoch epoch for a finite-horizon solutions.
#' @param ... further arguments are passed on.
#'
#' @returns `reward()` returns a vector of reward values, one for each belief if a matrix is specified.
#'
#' \item{state}{start state to calculate the reward for. if `NULL` then the start
#' state of model is used.}
#' @author Michael Hahsler
#' @examples
#' data("Maze")
#' Maze
#' gridworld_matrix(Maze)
#'
#' sol <- solve_MDP(Maze)
#' policy(sol)
#'
#' # reward for the start state s(3,1) specified in the model
#' reward(sol)
#'
#' # reward for starting next to the goal at s(1,3)
#' reward(sol, start = "s(1,3)")
#'
#' # expected reward when we start from a random state
#' reward(sol, start = "uniform")
#' @export
reward <- function(x, ...) {
  UseMethod("reward")
}

#' @rdname reward
#' @export
reward.MDP <- function(x,
                       start = NULL,
                       epoch = 1L,
                       ...) {
  is_solved_MDP(x, stop = TRUE)

  if (is.null(start)) {
    start <- start_vector(x)
  } else {
    start <- .translate_belief(start, x)
  }

  pol <- policy(x, epoch)

  if (is.null(pol$U)) {
    pol$U <- policy_evaluation(pol, x)
  }

  return(drop(crossprod(pol$U, start)))
}
