#' Solve MDPs using Linear Programming
#'
#' Solve discounted, infinite horizon MDPs via linear programming.
#'
#' @family solver
#' 
#' @details
#' A linear programming formulation was developed by Manne (1960) and further 
#' described by Puterman (1996). For the
#' optimal value function, the Bellman equation holds:
#'
#' \deqn{
#'  v^*(s) = \max_{a \in \mathcal{A}}\sum_{s' \in  \mathcal{S}} p(s, a, s') [ r(s, a, s') + \gamma v^*(s')]\; \forall a\in \mathcal{A}, s \in \mathcal{S}
#' }
#'
#' The maximization problem can reformulate as a minimization with 
#' a linear constraint for each state action pair. 
#' The optimal value function can be found by solving the 
#' following linear program:
#' \deqn{\text{min} \sum_{s\in S} v(s)}
#' subject to
#' \deqn{v(s) \ge \sum_{s' \in \mathcal{S}} p(s, a, s')[r(s, a, s') + \gamma v(s')],\; \forall a\in \mathcal{A}, s \in \mathcal{S}
#' }
#'
#' Note:
#' * The discounting factor has to be strictly less than 1.
#' * The used solver does not support infinity and a sufficiently large
#'   value needs to be used instead (see parameter `inf`).
#' * Additional parameters to to `solve_MDP` are passed on to [lpSolve::lp()].
#'
#' @references 
#' Manne, Alan. 1960. "On the Job-Shop Scheduling Problem." Operations Research 8 (2): 219-23. \doi{10.1287/opre.8.2.219}.
#' 
#' Puterman, Martin L. 1996. Markov decision processes: discrete stochastic dynamic programming. John Wiley & Sons.
#'
#' @examples
#' data(Maze)
#'
#' # we change the discount to 0.9 since LP is only implemented for discounted MDPs.  
#' maze_solved <- solve_MDP(Maze, discount = 0.9, method = "LP:LP", verbose = TRUE)
#' maze_solved
#' policy(maze_solved)
#' 
#' @inheritParams solve_MDP
#' @param method string; one of the following solution methods: `'LP'`
#' @param progress not supported by this solver.
#' @param horizon Only infinite-horizon MDPs with `horizon = Inf` are supported.
#' @param discount only undiscounted MDPs with `discount = 1` are supported.
#' @param inf value used for infinity when calling `lpSolve::lp()`. This
#'            should me much larger than the largest absolute
#'            reward in the model.
#' @param lpSolve_args a list with additional arguments passed on to `lpSolve::lp()`.
#' 
#' @inherit solve_MDP return
#' 
#' @export
solve_MDP_LP <- function(model,
                         method = "LP",
                         horizon = NULL,
                         discount = NULL,
                         inf = 1000,
                         lpSolve_args = list(),
                         ...,
                         matrix = NULL,
                         continue = FALSE,
                         verbose = FALSE,
                         progress = NULL
                         ) {
  .nodots(...)
  # currently ignored: matrix
  
  methods <- c("LP")
  method <- match.arg(method, methods)
  
  if (continue)
    stop("continue is not supported by LP methods!")
  
  model <- .prep_model(model, horizon, discount, matrix = FALSE, verbose, progress)
  
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
