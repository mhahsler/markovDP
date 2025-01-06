# Solve MDPs using Linear Programming

# For any given state, we have the assumption that the state's true value is:
# V^*(s) = r + \gamma \max_{a \in A}\sum_{s' \in S} P(s' | s,a) \cdot V^*(s')
#
# The non-linear max operation can be rewritten as a minimization with a
# linear inequality.
#
# min \sum_{s\in S} V(s)
# s.t.
# V(s) >=  R + \gamma\sum_{s' \in S} P(s' | s,a)*V(s'),\; \forall a\in A, s \in S
#
# find the smallest value of V(s) that matches this requirement,
# then that value would make exactly one of the constraints tight.


#' @rdname solve_MDP
#' @param inf value used for infinity when calling `lpSolve::lp()`. This
#'            should me much larger/smaller than the largest/smallest
#'            reward in the model.
#' @param lpSolve_args a list with additional arguments passed on to `lpSolve::lp()`.
#' @export
solve_MDP_LP <- function(model,
                         method = "lp",
                         horizon = NULL,
                         discount = NULL,
                         inf = 1000,
                         verbose = FALSE,
                         lpSolve_args = list(),
                         ...) {
  .nodots(...)
  # method is always "lp" and ignored
  
  ### horizon and discount rate
  if (!is.null(horizon)) {
    model$horizon <- horizon
  }
  if (is.null(model$horizon)) {
    model$horizon <- Inf
  }
  
  if (!is.null(discount)) {
    model$discount <- discount
  }
  if (is.null(model$discount)) {
    message("No discount factor specified. Using .9!")
    model$discount <- .9
  }
  
  if (is.finite(model$horizon)) {
    stop("method 'lp' can only be used for infinite horizon problems.")
  }
  
  gamma <- model$discount
  
  # TODO: For a better formulation of the undiscounted problem, see:
  # Lexing Ying and Yuhua Zhu, A note on optimization formulations of
  # Markov decision processes, Communications in Mathematical Sciences
  # 20(3), 727-745, 2022.
  # DOI: https://dx.doi.org/10.4310/CMS.2022.v20.n3.a5
  # https://arxiv.org/abs/2012.09417
  if (gamma >= 1) {
    warning("discount factor needs to be <1 for the used LP formulation. Using 0.999 instead of 1.")
    gamma <- 0.999
  }
  
  n_s <- length(S(model))
  n_a <- length(A(model))
  
  if (verbose) {
    cat("creating", n_s * n_a, "constraints ...")
  }
  
  # objective is sum_s V(s)
  obj <- rep(1, n_s)
  
  # enforce the Bellman equation for each state
  # V(s) \geq r + \gamma\sum_{s' \in S} P(s' | s,a)*V(s'),\; \forall a\in A, s \in S
  #
  # we can use the constraints V (I - gamma * P) >= PR for all a, s
  
  tm <- system.time({
    P <- do.call(rbind, transition_matrix(model, sparse = FALSE))
    
    # compute I - gamma * P for all a, s
    const_mat <- do.call(rbind, replicate(n_a, diag(n_s), simplify = FALSE)) - gamma * P
    
    R <- do.call(rbind, reward_matrix(model, sparse = FALSE))
    
    # fix R for lpSolve
    # lpSolve does not handle Inf
    R[R == +Inf] <- inf
    R[R == -Inf] <- -inf
    # lpSolve's simplex implementation requires all x > 0, but
    # state values can be negative! We add a second set of decision variables to
    # represent the negative values.
    neg_r <- any(R < 0)
    if (neg_r) {
      const_mat <- cbind(const_mat, -const_mat)
      obj <- c(rep(c(1, -1), each = n_s))
    }
    
    const_rhs <- rowSums(P * R)
  })
  
  if (verbose) {
    cat(" took", tm[1] + tm[2], "seconds.\n")
    cat("running LP solver ...")
  }
  
  tm <- system.time(solution <- do.call(lpSolve::lp, c(
    list(
      direction = "min",
      objective.in = obj,
      const.mat = const_mat,
      const.dir = rep(">=", length(const_rhs)),
      const.rhs = const_rhs
    ),
    lpSolve_args
  )))
  
  if (verbose) {
    cat(" took", tm[1] + tm[2], "seconds.\n")
  }
  
  if (solution$status) {
    stop("lp_solve did not solve the MDP successfully!")
  }
  
  U <- solution$solution
  
  # use the positive or the negative decision variable?
  if (neg_r) {
    U_neg <- U[-(1:n_s)]
    U <- U[1:n_s]
    U <- ifelse(U_neg < U, U, -U_neg)
  }
  
  pi <- greedy_policy(Q_values(model, U))
  
  model$solution <- list(
    method = "lp",
    policy = list(pi),
    converged = TRUE,
    solver_out = solution
  )
  
  model
}
