#' Regret of a Policy and Related Measures
#'
#' Calculates the regret and related measures for a policy relative to a benchmark policy.
#'
#' ## Regret
#' Regret for a policy \eqn{\pi} is defined as 
#' \deqn{v_\pi(s_0) - v_*(s_0),}
#' where \eqn{v_\pi(s_0)} represents the expected long-term
#' state value for following policy \eqn{\pi} and the starting
#' in state \eqn{s_0} (or a start distribution).
#' The relative regret is calculated as 
#' \deqn{\frac{v_\pi(s_0) - v_*(s_0)}{v_*(s_0)}.}
#' 
#' Note that for regret, usually the optimal policy \eqn{\pi^*} is used as the benchmark.
#' Since the optimal policy may not be known, regret relative to the best known policy can be used.
#'
#' Regret is only valid with converged value functions. This means that either the
#' solver has converged, or the value function was estimated for the policy using
#' converged [policy_evaluation()].
#'
#' ## Action Discrepancy
#' The action discrepancy 
#' measures the difference between two policies as the number 
#' of states for which the prescribed action in the policies differs. 
#' Often, a policy is compared to the best known policy called the benchmark
#' policy.
#' 
#' Some times two actions are equivalent (have the same q-value) and 
#' the algorithm breaks the tie randomly. The implementation accounts for this 
#' case.
#' 
#' The action discrepancy can be calculated as a proportion of different actions
#' or be weighted by the state visit probability given the benchmark policy. Both 
#' weighted and proportional action discrepancy is scaled in \eqn{[0, 1]}.
#'
#' ## Root Mean Squared Value Error
#' The root mean value error 
#' \deqn{\sqrt{\text{VE}} = \sqrt{||v_\pi - v_*||^2}} 
#' is the sum of the squared differences of
#' state values between a solution's value function and the optimal value function. 
#' For \eqn{v_*}, the
#' value function of the benchmark solution is used.
#' Related measures like MSVE (means squared value error), 
#' MAVE (means absolute value error), 
#' AVE (absolute value error) and VE (value error)
#' are also provided.
#' 
#' The error can also be weighted by the state visit probability 
#' given the benchmark policy. This results in the expected 
#' error with respect to the state visit distribution of the benchmark policy. 
#' This may be important to evaluate methods that focus 
#' only on estimating the value function for states that are actually visited using the 
#' policy. 
#'
#' @family MDP
#' @family policy
#'
#' @param policy a solved MDP containing the policy to calculate the regret for.
#' @param benchmark a solved MDP with the (optimal) policy. Regret is calculated relative to this
#'    policy.
#' @param start start state distribution. If NULL then the start state of the `benchmark` is used.
#' @param relative logical; should the relative regret (regret divided by the reward of the benchmark)
#' be calculated?
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
#'
#' # a manual policy (go up and in some squares to the right)
#' acts <- rep("up", times = length(Maze$states))
#' names(acts) <- Maze$states
#' acts[c("s(1,1)", "s(1,2)", "s(1,3)")] <- "right"
#'
#' sol_manual <- add_policy(Maze, manual_policy(Maze, acts, estimate_V = TRUE))
#' 
#' # compare the policies side-by-side
#' cbind(opt = policy(sol_optimal), manual = policy(sol_manual))
#'
#' # the regret is very small. It is about 4.8% of the optimal reward
#' regret(sol_manual, benchmark = sol_optimal)
#' regret(sol_manual, benchmark = sol_optimal, relative = TRUE)
#'
#' # The number of different actions (excluding equivalent actions) is 3.
#' # This about 27% of the actions in the policy. 
#' action_discrepancy(sol_manual, benchmark = sol_optimal)
#' action_discrepancy(sol_manual, benchmark = sol_optimal, proportion = TRUE)
#' 
#' # Weighted by the probability that a state will be visited shows that
#' only 2.3% of the time a different action would be used.
#' action_discrepancy(sol_manual, benchmark = sol_optimal, weighted = TRUE)
#'
#' value_error(sol_manual, benchmark = sol_optimal, type = "VE")
#' value_error(sol_manual, benchmark = sol_optimal, type = "MAVE")
#' value_error(sol_manual, benchmark = sol_optimal, type = "RMSVE")
#' 
#' # Weighting shows that the expected MAVE (expectation taken over the 
#' # benchmark policy) is rather small.
#' value_error(sol_manual, benchmark = sol_optimal, type = "MAVE", weighted = TRUE)
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
                       relative = FALSE,
                       ...) {
  rb <- reward(benchmark, start, ...)
  rp <- reward(policy, start, ...)
  
  if (relative)
    (rb - rp) / rb
  else
    rb - rp
}

#' @rdname regret
#' @param weighted logical; should mismatched actions or state value errors 
#'   be weighted by the state visit probability for the benchmark? Rarely or never
#'   visited states will have now less influence on the measure.
#' @param proportion logical; should the action discrepancy be reported as a proportion
#'    of states with a different action.
#' @export
action_discrepancy <- function(policy, benchmark, weighted = FALSE, proportion = FALSE) {
  if (weighted && proportion)
    stop("You can only use weighted or proportion!")
  
  is_solved_MDP(benchmark)
  is_solved_MDP(policy)
  
  pi <- policy(policy)$action
  Q <- Q_values(benchmark)
  
  if (weighted)
    weight <- visit_probability(benchmark)
  else 
    weight <- 1
  
  # account for ties in the Q matrix
  discrepancy <- sum(weight * (Q[cbind(seq_len(nrow(Q)), pi)] != apply(Q, MARGIN = 1, max)))
  
  if (proportion)
    discrepancy <- discrepancy / length(pi)
  
  discrepancy
}

#' @rdname regret
#' @param type type of error root mean square value error (`"RMSVE"`), mean square value error (`"MSVE"`), 
#'   mean absolute value error (`"MAVE"`), absolute value error vector (`"AVE"`),  value error vector (`"VE"`).
#' @export
value_error <- function(policy, benchmark, type = "RMSVE", weighted = FALSE) {
  is_solved_MDP(benchmark)
  is_solved_MDP(policy)
  
  type <- match.arg(type, c("RMSVE", "MSVE", "MAVE", "AVE", "VE"))
  
  ve <- abs(policy(policy)$V - policy(benchmark)$V)
  
  if (weighted)
    ve <- ve * visit_probability(benchmark)
  
  switch(type,
         VE = ve,
         AVE= abs(ve),
         MAVE = mean(abs(ve)),
         MSVE = mean(ve^2),
         RMSVE = sqrt(mean(ve^2))
  )
}
