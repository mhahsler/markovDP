#' Calculate the Expected Reward of a Policy
#'
#' This function calculates the expected total reward for an MDP policy
#' given a start state (distribution). The value is calculated using the value
#' function stored in the MDP solution.
#'
#' The reward is typically calculated using the value function
#' of the solution. If these are not available, then [sample_MDP()] is
#' used instead with a warning.
#'
#' @family policy
#'
#' @param model a solved [MDP] object.
#' @param start specification of the current state (see argument start
#' in [MDP] for details). By default the start state defined in
#' the model as start is used. Multiple states can be specified as rows in a matrix.
#' @param method `"solution"` uses the converged value function stored in the solved model, 
#'               `"policy_evaluation"` estimates the value function, and `"sample"`
#'               calculates the average reward by sampling episodes from the model.
#' @param ... further arguments are passed on to [policy_evaluation()] or [sample_MDP()].
#'
#' @returns `reward()` returns a vector of reward values, one for each belief if a matrix is specified.
#'
#' \item{state}{start state to calculate the reward for. if `NULL` then the start
#' state of model is used.}
#' @author Michael Hahsler
#' @examples
#' data("Maze")
#' Maze
#' gw_matrix(Maze)
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
#' # expected reward when we start from a random state as returned from the solver
#' reward(sol, start = "uniform")
#' 
#' # estimate the reward using sampling following the policy
#' reward(sol, method = "sample", start = "uniform", n = 10000, horizon = 1000)
#' @export
reward <- function(model, ...) {
  UseMethod("reward")
}

#' @rdname reward
#' @export
reward.MDP <- function(model,
                       start = NULL,
                       method = "solution",
                       ...) {
  method <- match.arg(method, c("solution", "policy_evaluation", "sample"))
  start <- start_vector(model, start = start)
 
  if (method == "solution" && !is_converged_MDP(model)) {
    method <- "policy_evaluation"
    warning("model does not contain a converged solution. Using policy evaluation to obtain the value function.")
  }
   
  if (method == "solution") {
    r <- sum(policy(model)$V * start)
  }
  
  else if (method == "policy_evaluation") {
    r <- sum(policy_evaluation(model, policy(model), ...) * start)
  }
  
  else if (method == "sample") {
    r <- sample_MDP(model, start = start, ...)$avg_reward
  }
  
  else
    stop("Unknown method!")
  
  r
}
