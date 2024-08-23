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
#' @family MDP
#'
#' @param policy a solved MDP containing the policy to calculate the regret for.
#' @param benchmark a solved MDP with the (optimal) policy. Regret is calculated relative to this
#'    policy.
#' @param start start state distribution. If NULL then the start state of the `benchmark` is used.
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
#' sol_manual <- add_policy(Maze, manual_policy(Maze, acts))
#' policy(sol_manual)
#'
#' regret(sol_manual, benchmark = sol_optimal)
#' @export
regret <- function(policy, benchmark, start = NULL) {
  UseMethod("regret")
}

#' @export
regret.MDP <- function(policy, benchmark, start = NULL) {
  if (!inherits(benchmark, "MDP") || !is_solved_MDP(benchmark)) {
    stop("benchmark needs to be a solved MDP.")
  }

  if (!inherits(policy, "MDP") || !is_solved_MDP(policy)) {
    stop("policy needs to be a solved MDP.")
  }

  if (is.null(start)) {
    start <- which(start_vector(benchmark) == 1)
  }

  if (is.character(start)) {
    start <- which(benchmark$states == start)
  }

  if (length(start) != 1L) {
    stop("A single start state needs to be specified!")
  }

  r_bench <- policy(benchmark)$U[start]
  r_pol <- policy(policy)$U[start]

  r_bench - r_pol
}
