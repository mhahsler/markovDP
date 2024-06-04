# Solve MDPs using Linear Programming

# For any given state, we have the assumption that the state's true value is:
# V^*(s) = r + \gamma \max_{a \in A}\sum_{s' \in S} T(s' | s,a) \cdot V^*(s')
#
#
# min \sum_{s\in S} V(s)
# s.t.
# V(s) >=  R + \gamma\sum_{s' \in S} P(s' | s,a)*V(s'),\; \forall a\in A, s \in S
#
# find the smallest value of V(s) that matches this requirement,
# then that value would make exactly one of the constraints tight.


#' @rdname solve_MDP
#' @export
solve_MDP_LP <- function(model, method = NULL, verbose = FALSE, ...) {
  # method is always "lp" and ignored

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
    warning("discount factor needs to be <1 for LP. Using 0.999 instead of 1.")
    gamma <- 0.999
  }

  n_s <- length(model$states)
  n_a <- length(model$actions)

  # objective is sum_s V(s)
  obj <- rep(1, n_s)

  # enforce the Bellman equation for each state
  # V(s) \geq r + \gamma\sum_{s' \in S} P(s' | s,a)*V(s'),\; \forall a\in A, s \in S
  #
  # we can use the constraints V (I - GAMMA * T) >= TR for all a, s

  # T(s,a, s'): maps actions + start states -> end states
  T <- matrix(NA_real_, nrow = n_a * n_s, ncol = n_s)
  i <- 1L
  for (a in seq(n_a)) {
    for (s in seq(n_s)) {
      T[i, ] <- transition_matrix(model, action = a, start.state = s)
      i <- i + 1L
    }
  }
  rownames(T) <- paste(rep(model$actions, each = n_s), rep(model$states, n_a))
  colnames(T) <- model$states
  
  const_mat <- do.call(rbind, replicate(n_a, diag(n_s), simplify = FALSE)) - gamma * T
  const_dir <- rep(">=", n_a * n_s)

  # make reward compatible into a S x S matrix
  R <- sapply(model$actions, FUN = function(a) {
    t(do.call(cbind, reward_matrix(model, action = a)))
  }, simplify = FALSE)

  # lpSolve's simplex implementation requires all x > 0, but 
  # state values can be negative! We add a second set of decision variables to 
  # represent the negative values.
  neg_r <- min(sapply(R, min)) < 0
  if (neg_r) {
    const_mat <- cbind(const_mat , -const_mat)
    obj <- c(rep(c(1, -1), each = n_s))
  }

  # sum_sp (T * R): maps actions + start states -> expected reward
  TR <- sapply(model$actions, FUN = function(a) {
    rowSums(transition_matrix(model, action = a) * R[[a]])
  })
  TR <- as.vector(TR)
  names(TR) <- paste(rep(model$actions, each = n_s), rep(model$states, n_a))
  const_rhs <- TR

  if (verbose) {
    cat("running LP solver ...")
  }

  tm <- system.time(solution <- lpSolve::lp(
    direction = "min",
    objective.in = obj,
    const.mat = const_mat,
    const.dir = const_dir,
    const.rhs = const_rhs,
    ...
  ))

  if (verbose) {
    cat(" took", tm[1] + tm[2], "seconds.\n")
  }

  if (solution$status) {
    stop("lp_solve did not solve the MDP successfully!")
  }

  U <- solution$solution
  
  # use the positive or the negative decision variable.
  if(neg_r) {
    U_neg <- U[-(1:n_s)]
    U <- U[1:n_s]
    neg <- U_neg > U
    U[neg] <- -U_neg[neg]
  }
  
  
  pi <- greedy_policy(q_values(model, U))

  model$solution <- list(
    method = "lp",
    policy = list(pi),
    converged = NA,
    solver_out = solution
  )

  model
}
