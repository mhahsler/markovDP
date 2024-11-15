#' State Visit Probability 
#'
#' Calculates the state visit probability (the modified stationary distribution) 
#' when following a policy from the start state.
#'
#' The visit probability is the stationary distribution for the transition matrix
#' induced by the policy. To account for absorbing states, we modify the
#' transition matrix by setting all outgoing probabilities from absorbing states 
#' to 0. 
#' 
#' ## Power iteration
#' 
#' The stationary distribution can be estimated as the sum of multiplying the 
#' start distribution repeatedly with the modified transition matrix 
#' induced by the policy. We stop multiplying when 
#' the largest difference between entries in the two consecutive vectors is 
#' less then the extra parameter:
#' 
#' *  `min_err` stop criterion for the power iteration (default: `1e-6`).  
#' 
#' The resulting vector is normalized to probabilities.
#'
#' ## Sample method
#' The stationary distribution is calculated using `n` random walks. The state visit 
#' counts are normalized to a probabilities.
#' Additional parameters are:
#' 
#' * `n` number of random walks (default `1000`).
#' * `horizon` maximal horizon used to stop a random walk if it has 
#'    not reached an absorbing state.
#'
#' @family policy
#'
#' @param model a solved [MDP] object.
#' @param start specification of the start distribution. If missing the specification in
#'    `model` is used.
#' @param pi the used policy. If missing the policy in
#'    `model` is used.
#' @param method calculate the modified stationary distribution using `"power"` 
#'     (power iteration) or `"sample"` (trajectory sampling).  
#' @param min_err repeats multiplying till the largest difference between 
#'     two consecutive vectors is less then `min_err`.
#'
#' @returns a visit probability vector over all states.
#'
#' @author Michael Hahsler
#' @examples
#' data("Maze")
#' Maze
#' 
#' sol <- solve_MDP(Maze)
#' visit_probability(sol)
#' 
#' # gw_matrix also can calculate the visit_probability.
#' gw_matrix(sol, what = "visit_probability")
#'
# @export
visit_probability <- function(model, pi = NULL, start = NULL, method = "power", ...) {
    method <- match.arg(method, c("matrix", "sample")) 
  
    switch(method,
      matrix = visit_probability_power(model, pi, start, ...),
      sample = visit_probability_sample(model, pi, start, ...),
    )
  }

visit_probability_power <- function(model, pi = NULL, start = NULL, max_err = 1e-6) {
  if (is.null(pi))
    pi <- policy(model)
  
  start <- start_vector(model, start = start, sparse = FALSE)
  
  P <- induced_transition_matrix(model, pi)
  absorbing <- absorbing_states(model, sparse = "index")
  
  P[absorbing, ] <- 0
  
  X <- start %*% P
  S <- start + X

  err <- Inf
  while (err > max_err) { 
    Xp <- X %*% P
    err <- max(abs(X-Xp))
    X <- Xp
    S <- S + X 
  }
  
  drop(S/sum(S))
}

visit_probability_sample <- function(model, pi = NULL, start = NULL, n = 1000, horizon = NULL) {
  if (!is.null(pi))
    model <- add_policy(model, pi)
  
  start <- start_vector(model, start = start, sparse = FALSE)
  
  samp <- sample_MDP(model, n, start = start, horizon = horizon, epsilon = 0)
  samp$state_cnt / sum(samp$state_cnt)
}
