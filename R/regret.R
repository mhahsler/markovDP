#' Regret of a Policy and Related Measures
#'
#' Calculates the regret and related measures for a policy relative to a benchmark policy.
#'
#' ## Regret 
#' Regret is defined as \eqn{V^{\pi^*}(s_0) - V^{\pi}(s_0)} with \eqn{V^\pi} representing the expected long-term
#' state value (represented by the value function) given the policy \eqn{\pi} and the start
#' state \eqn{s_0}.
#'
#' Note that for regret, usually the optimal policy \eqn{\pi^*} is used as the benchmark.
#' Since the optimal policy may not be known, regret relative to the best known policy can be used.
#' 
#' Regret is only valid with converged value functions. This means that either the
#' solver has converged, or the value function was estimated for the policy using
#' converged [policy_evaluation()].
#' 
#' ## Action Discrepancy
#' The action discrepancy is the number of states in the policy for which the 
#' prescribed action in the policy differs. This implementation calculates
#' the Q matrix for the benchmark to make sure that actions with the same Q-value
#' are considers as correct. 
#' 
#' ## Mean (squared) Value Error 
#' The mean value error is the sum of the squared differences in 
#' state values between the two solutions. 
#'
#' @family MDP
#'
#' @param policy a solved MDP containing the policy to calculate the regret for.
#' @param benchmark a solved MDP with the (optimal) policy. Regret is calculated relative to this
#'    policy.
#' @param start start state distribution. If NULL then the start state of the `benchmark` is used.
#' @param  policy_eval logical; run policy evaluation to re-estimate state values.
#' @param proportion logical; should the action discrepancy be reported as a proportion
#'    of states with a different action.
#' @param square logical; should the value error be squared?
#' @param ... further arguments are passed on to [policy_evaluation()].
#' 
#' @returns 
#'    * `regret()` returns the regret as a difference of expected long-term rewards.
#'    * `action_discrepancy()` returns the number or proportion of diverging actions.
#'    * `mean_value_error()` returns the mean squared or absolute difference in the
#'        value function.
#'
#' @author Michael Hahsler
#' @examples
#' data(Maze)
#'
#' sol_optimal <- solve_MDP(Maze)
#' policy(sol_optimal)
#'
#' # a manual policy (go up and in some squares to the right)
#' acts <- rep("up", times = length(Maze$states))
#' names(acts) <- Maze$states
#' acts[c("s(1,1)", "s(1,2)", "s(1,3)")] <- "right"
#' 
#' sol_manual <- add_policy(Maze, manual_policy(Maze, acts, estimate_U = TRUE))
#' policy(sol_manual)
#'
#' regret(sol_manual, benchmark = sol_optimal)
#' 
#' action_discrepancy(sol_manual, benchmark = sol_optimal)
#' 
#' mean_value_error(sol_manual, benchmark = sol_optimal)
#' @export
regret <- function(policy, benchmark, start = NULL, run_policy_eval = TRUE, ...) {
  UseMethod("regret")
}

#' @export
regret.MDP <- function(policy, benchmark, start = NULL, policy_eval = TRUE, ...) {
  is_solved_MDP(benchmark)
  is_solved_MDP(policy)
  
  if (!policy_eval &&
      (!is_converged_MDP(benchmark) || !is_converged_MDP(policy)))
    warning("The MDPs need to have converged if no policy evaluation is used.")

  start <- start_vector(benchmark, start = start)
  
  if (policy_eval) {
    U_bench <- policy_evaluation(benchmark, policy(benchmark), ...)
    U_pol <- policy_evaluation(benchmark, policy(policy), ...)
  } else {
    U_bench <- policy(benchmark)$U
    U_pol <- policy(policy)$U
  }
  
  
  r_bench <- sum(U_bench * start)
  r_pol <- sum(U_pol * start)

  r_bench - r_pol
}

#' @rdname regret
#' @export
action_discrepancy <- function(policy, benchmark, proportion = FALSE) {
  is_solved_MDP(benchmark)
  is_solved_MDP(policy)
  
  pi <- policy(policy)$action
  Q <- q_values(benchmark)
  
  # account for ties in the Q matrix
  discrepancy <- sum(Q[cbind(seq_len(nrow(Q)), pi)] != apply(Q, MARGIN = 1, max))
  
  if (proportion)
    discrepancy <- discrepancy / length(pi)
  
  discrepancy
}

#' @rdname regret
#' @export
mean_value_error <- function(policy, benchmark, square = TRUE) {
  is_solved_MDP(benchmark)
  is_solved_MDP(policy)
  
  ve <- abs(policy(policy)$U - policy(bench)$U)
  if (square)
    ve <- ve^2
  mean(ve)
}
