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
solve_MDP_LP <- function(model, method = NULL, ...) {
  # method is always "lp" and ignored
  
  if (is.finite(model$horizon))
    stop("mtehod 'lp' and only be used for infinite horizon problems.")
  
  gamma <- model$discount
  if (gamma >= 1) {
    warning("discount factor needs to be <1 for LP. Using 0.999.")
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
  T <- NULL
  for(a in model$actions)
    for(s in model$states)
      T <- rbind(T, transition_matrix(model, action = a, start.state = s))
  rownames(T) <- paste(rep(model$actions, each = n_s), rep(model$states, n_a))

  const_mat <- do.call(rbind, replicate(n_a, diag(n_s), simplify=FALSE)) - gamma * T
  const_dir <- rep(">=", n_a * n_s)

  # make reward compatible into a S x S matrix
  R <- sapply(model$actions, FUN = function(a) 
    t(do.call(cbind, reward_matrix(model, action = a))), simplify = FALSE)

  # lpSolve requires all x > 0! 
  # we make the reward non-negative
  R_shift <- FALSE
  R_min <- min(sapply(R, min))
  if (R_min < 0) {
    R_shift <- TRUE
    warning("negative rewards. Shifting rewards for lpSolve.")
    for (i in seq_along(R))
      R[[i]] <- R[[i]] - R_min
  }
    
  # sum_sp (T * R): maps actions + start states -> expected reward
  TR <- sapply(model$actions, FUN = function(a) 
    rowSums(transition_matrix(model, action = a) * R[[a]]))
  TR <- as.vector(TR)
  names(TR) <- paste(rep(model$actions, each = n_s), rep(model$states, n_a))
  const_rhs <- TR

  solution <- lpSolve::lp(direction = "min", 
                 objective.in = obj, 
                 const.mat = const_mat, 
                 const.dir = const_dir, 
                 const.rhs = const_rhs,
                 ...)
  
  if (solution$status)
    stop("lp_solve did not solve the MDP successfully!")
  
  U <- solution$solution
  pi <- greedy_MDP_policy(q_values_MDP(model, U))
 
  # update U with unshifted values
  if (R_shift) 
    pi$U <- MDP_policy_evaluation(pi, model,verbose = TRUE)
   
  model$solution <- list(
    method = "lp",
    policy = list(pi),
    converged = NA,
    solver_out = solution
  )
  
  model
}

