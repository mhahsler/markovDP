# Policy evaluation using Linear Programming

# For any given state, we have the assumption that the state's true value
#   given policy pi is:
# V^pi(s) = r + \gamma \sum_{s' \in S} p(s' | s, pi(s)) \cdot V^pi(s')
#
# min \sum_{s\in S} V(s)
# s.t.
# V(s) >=  r + \gamma\sum_{s' \in S} p(s' | s, pi(s))*V(s'),\; \forall s \in S


#' @rdname policy_evaluation
#' @details The LP implementation only works for infinite-horizon problems with 
#' a discount factor <1. It solves the following LP:
#' 
#' \deqn{\min \sum_{s\in S} v(s)}
#' s.t.
#' \deqn{v(s) \ge \sum_{s' \in S} p(s' | s, \pi(s)) [r(s, \pi(s), s') + \gamma v(s')],\; \forall s \in S         
#'  }        
#'          
#' @param inf value used to replace infinity for `lpSolve::lp()`.
#' @param ... further arguments are ignored
#' @export
policy_evaluation_LP <- function(model,
                                 pi = NULL,
                                 inf = 1000, 
                                 ...,
                                 progress = FALSE,
                                 verbose = FALSE) {
 
  # progress not supported! 
  
  if (is.finite(model$horizon)) {
    stop("method 'lp' can only be used for infinite horizon problems.")
  }
  
  gamma <- model$discount %||% 1
  
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
    cat("creating constraints ...\n")
  }
  
  # to enforce the Bellman equation for each state
  # V(s) \geq r + \gamma\sum_{s' \in S} P(s' | s,a)*V(s'),\; \forall  s \in S
  #
  # we can use the constraints V (I - gamma * T) >= TR for all s
  
  # objective is sum_s V(s)
  obj <- rep(1, n_s)
  
  P <- induced_transition_matrix(model, pi)
  const_mat <- diag(n_s) - gamma * P
  
  R <- induced_reward_matrix(model, pi)
  # lsSolve does not handle Inf
  R[R == +Inf] <- inf
  R[R == -Inf] <- -inf
  
  # lpSolve's simplex implementation requires all x > 0, but
  # state values can be negative! We add a second set of decision variables to
  # represent the negative values.
  neg_r <- any(R < 0)
  if (neg_r) {
    const_mat <- cbind(const_mat, -const_mat)
    obj <- c(obj, -obj)
  }
  
  PR <- rowSums(P * R)
  
  if (verbose) {
    cat("running LP solver ...")
  }
  
  tm <- system.time(
    solution <- lpSolve::lp(
      direction = "min",
      objective.in = obj,
      const.mat = const_mat,
      const.dir = rep(">=", length(PR)),
      const.rhs = PR,
      ...
    )
  )
  
  if (verbose) {
    cat(" took", tm[1] + tm[2], "seconds.\n")
  }
  
  if (solution$status) {
    stop("lp_solve did not solve the MDP successfully!")
  }
  
  U <- solution$solution
  
  # use the positive or the negative decision variable.
  if (neg_r) {
    U_neg <- U[-(1:n_s)]
    U <- U[1:n_s]
    U <- ifelse(U_neg < U, U, -U_neg)
  }
  
  U
}
