#' Calculate the Regret of a Policy
#'
#' Calculates the regret of a policy relative to a benchmark policy.
#'
#' Regret is defined as \eqn{V^{\pi^*}(s_0) - V^{\pi}(s_0)} with \eqn{V^\pi} representing the expected long-term
#' state value (represented by the value function) given the policy \eqn{\pi} and the start
#' state \eqn{s_0}.
#'
#' Note that for regret usually the optimal policy \eqn{\pi^*} is used as the benchmark.
#' Since the optimal policy may not be known, regret relative to the best known policy can be used.
#' 
#' If the policy has 
#'
#' @family MDP
#'
#' @param policy a solved MDP containing the policy to calculate the regret for.
#' @param benchmark a solved MDP with the (optimal) policy. Regret is calculated relative to this
#'    policy.
#' @param start start state distribution. If NULL then the start state of the `benchmark` is used.
#' @param  run_policy_eval logical; run policy evaluation to re-estimate state values.
#' @param ... further arguments are passed on to [policy_evaluation()].
#' 
#' @return the regret as a difference of expected long-term rewards.
#'
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
#' U <- policy_evaluation(Maze, manual_policy(Maze, acts))
#' 
#' sol_manual <- add_policy(Maze, manual_policy(Maze, acts, U = TRUE))
#' policy(sol_manual)
#'
#' regret(sol_manual, benchmark = sol_optimal)
#' @export
regret <- function(policy, benchmark, start = NULL, run_policy_eval = TRUE, ...) {
  UseMethod("regret")
}

#' @export
regret.MDP <- function(policy, benchmark, start = NULL, run_policy_eval = TRUE, ...) {
  if (!inherits(benchmark, "MDP") || !is_solved_MDP(benchmark)) {
    stop("benchmark needs to be a solved MDP.")
  }

  if (!inherits(policy, "MDP") || !is_solved_MDP(policy)) {
    stop("policy needs to be a solved MDP.")
  }

  start <- start_vector(benchmark, start = start)
  
  if (run_policy_eval) {
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
