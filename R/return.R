#' Calculate the Expected Return of a Policy
#'
#' This function calculates the expected total return for an MDP policy
#' given a start state (distribution). The value is calculated using the value
#' function stored in the MDP solution.
#'
#' The return is typically calculated using the value function
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
#'               calculates the average return by sampling episodes from the model.
#' @param ... further arguments are passed on to [policy_evaluation()] or [sample_MDP()].
#'
#' @returns `expected_return()` returns a vector of returns, one for each start state if a matrix is specified.
#'
#' \item{state}{start state to calculate the return for. If `NULL` then the start
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
#' # return for the start state s(3,1) specified in the model
#' expected_return(sol)
#'
#' # return for starting next to the goal at s(1,3)
#' expected_return(sol, start = "s(1,3)")
#'
#' # expected return when we start from a random state as returned from the solver
#' expected_return(sol, start = "uniform")
#' 
#' # estimate the return using sampling following the policy
#' expected_return(sol, method = "sample", start = "uniform", n = 10000, horizon = 1000)
#' @export
expected_return <- function(model, ...) {
  UseMethod("expected_return")
}

#' @rdname expected_return
#' @export
expected_return.MDP <- function(model,
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
    r <- sample_MDP(model, start = start, ...)$avg_return
  }
  
  else
    stop("Unknown method!")
  
  r
}
