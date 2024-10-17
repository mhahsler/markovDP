#' Regret of a Policy and Related Measures
#'
#' Calculates the regret and related measures for a policy relative to a benchmark policy.
#'
#' ## Regret
#' Regret for a policy \eqn{\pi} is defined as 
#' \eqn{V_\pi(s_0) - v_*(s_0)}.
#' \eqn{V_\pi(s_0)} representing the expected long-term
#' state value for following policy \eqn{\pi} and the starting
#' in state \eqn{s_0} (or a start distribution).
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
#' the Q matrix for the benchmark to make sure that actions with the same Q value
#' are considers as correct.
#'
#' ## Root Mean Squared Value Error
#' The root mean value error \eqn{\sqrt{\text{VE}} = \sqrt{||V - v_*||^2}} 
#' is the sum of the squared differences in
#' state values between the two solution and the optimal value function. For \eqn{v_*}, the
#' value function of the benchmark solution is used.
#'
#' @family MDP
#'
#' @param policy a solved MDP containing the policy to calculate the regret for.
#' @param benchmark a solved MDP with the (optimal) policy. Regret is calculated relative to this
#'    policy.
#' @param start start state distribution. If NULL then the start state of the `benchmark` is used.
#' @param ... further arguments are passed on to [reward()].
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
#' sol_manual <- add_policy(Maze, manual_policy(Maze, acts, estimate_V = TRUE))
#' policy(sol_manual)
#'
#' regret(sol_manual, benchmark = sol_optimal)
#'
#' action_discrepancy(sol_manual, benchmark = sol_optimal)
#'
#' value_error(sol_manual, benchmark = sol_optimal, type = "VE")
#' value_error(sol_manual, benchmark = sol_optimal, type = "MAVE")
#' value_error(sol_manual, benchmark = sol_optimal, type = "RMSVE")
#' @export
regret <- function(policy,
                   benchmark,
                   start = NULL,
                   ...) {
  UseMethod("regret")
}

#' @export
regret.MDP <- function(policy,
                       benchmark,
                       start = NULL,
                       ...) {
   reward(benchmark, start, ...) - reward(policy, start, ...)
}

#' @rdname regret
#' @param proportion logical; should the action discrepancy be reported as a proportion
#'    of states with a different action.
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
#' @param type type of error root mean square value error (`"RMSVE"`), mean square value error (`"MSVE"`), 
#'   mean absolute value error (`"MAVE"`), absolute value error vector (`"AVE"`),  value error vector (`"VE"`).
#' @export
value_error <- function(policy, benchmark, type = "RMSVE") {
  is_solved_MDP(benchmark)
  is_solved_MDP(policy)
  
  type <- match.arg(type, c("RMSVE", "MSVE", "MAVE", "AVE", "VE"))
  
  ve <- abs(policy(policy)$V - policy(benchmark)$V)
  
  switch(type,
         VE = ve,
         AVE= abs(ve),
         MAVE = mean(abs(ve)),
         MSVE = mean(ve^2),
         RMSVE = sqrt(mean(ve^2))
  )
}
