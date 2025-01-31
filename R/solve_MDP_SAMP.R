#' Solve MDPs using Random-Sampling
#'
#' Solve MDPs via random-sampling. Implemented is Sutton and Barto's one-step 
#' random-sampling one-step tabular Q-planning.
#'
#' @family solver
#' 
#' @details
#' Random-sample one-step tabular Q-planning is a simple, not very effective, 
#' planning method shown as an illustration in Chapter 8 of
#' Sutton and Barto (2018). It randomly selects a
#' state/action pair and samples the following state \eqn{s'} to 
#' perform a single a one-step update:
#' 
#' \deqn{Q(s,a) \leftarrow Q(s, a) + \alpha * (R + \gamma \max_a'(Q(s', a')) - Q(s, a))}
#' 
#' @references
#' 
#' Sutton, Richard S., and Andrew G. Barto. 2018. Reinforcement Learning: An Introduction. Second. The MIT Press. [http://incompleteideas.net/book/the-book-2nd.html](http://incompleteideas.net/book/the-book-2nd.html).
#' 
#' @inheritParams solve_MDP
#' @param method string; one of the following solution methods: `'q_planning'`
#' @param alpha step size as a function of the time step `t` and the number of times
#'   the respective Q-value was updated `n` or a scalar. For expected Sarsa, alpha is
#'   often set to 1.
#' @param n number of updates performed.
#' @param Q an initial state-action value matrix. By default an all 0 matrix is 
#'        used.
#' @inherit solve_MDP return 
#' 
#' @export
solve_MDP_SAMP <-
  function(model,
           method = "q_planning",
           horizon = NULL,
           discount = NULL,
           alpha = function(t, n) min(10/n, 1),
           n = 1000,
           Q = NULL,
           ...,
           matrix = TRUE,
           continue = FALSE,
           progress = TRUE,
           verbose = FALSE) {
    .nodots(...)
    PROGRESS_INTERVAL <- 100
    
    method <-
      match.arg(method, c("q_planning"))
    
    if (!is.null(horizon))
      warning("q_planning does not use horizon. The specified horizon is ignored.")
    
    model <- .prep_model(model, horizon, discount, matrix, verbose, progress)
    
    if (!is.function(alpha))
      alpha_val <- alpha
    
    # Initialize Q
    if (continue) {
      if (is.null(model$solution$Q) || is.null(model$solution$Q_N))
        stop("model solution does not contain a Q matrix or the N count matrix to continue from!")
      Q <- model$solution$Q
      Q_N <- model$solution$Q_N
    } else if (is.null(Q)) {
      Q <- Q_zero(model)
      Q_N <-
        matrix(0L,
               nrow = nrow(Q),
               ncol = ncol(Q),
               dimnames = dimnames(Q)
        )
    }
    
    # return unconverged result when interrupted
    on.exit({
      if (progress) {
        pb$tick(0)
        pb$terminate()
      }
      
      if (t < n)
        warning("Manual interupt: MDP solver stopped at try ", t)
      
      if (verbose) {
        cat("\nTerminated after try:", t, "\n")
      }
      
      model$solution <- list(
        method = method,
        alpha = alpha,
        n = n,
        Q = Q,
        Q_N = Q_N,
        converged = NA,
        policy = list(greedy_policy(Q))
      )
      return(model)
    })
    
    if (progress)
      pb <- my_progress_bar(ceiling(n/PROGRESS_INTERVAL) + 1, name = "solve_MDP")
    
    if (verbose) {
      cat("Running q_planning (sampling)")
      cat("\nalpha:       ", deparse(alpha))
      cat("\nn:           ", n)
      
      cat("\nInitial Q (first 20 max):\n")
      print(head(Q, n = 20))
      
      cat("\nInitial Q_N (first 20 max):\n")
      print(head(Q_N, n = 20))
      cat("\n")
    }
    
    S <- S(model)
    A <- A(model)
    start <- start_vector(model, sparse = FALSE)
    gamma <- model$discount
    
    # loop through tries
    t <- 0L
    while (t < n) {
      t <- t + 1L
      if (!(t %% PROGRESS_INTERVAL) && progress)
        pb$tick()
      
      # sample a state/action pair
      s <- sample(S, 1L)
      a <- sample(A, 1L)
      
      sp_r <- act(model, s, a, fast = TRUE)
      sp <- sp_r$state_prime
      r <- sp_r$r
      
      # NOTE: we use as the default alpha = 1/N(s,a)
      #       then the expected error for each Q value is sigma_TD/sqrt(N(s,a))
      #           sigma_TD is the sd of R+gamma max_a Q(s', a)
     
      
      # update Q and Q_N
      Q_N[s, a] <- Q_N[s, a] + 1L
      
      if (is.function(alpha))
        alpha_val <- alpha(t, Q_N[s, a])
      
      if (verbose > 1) 
        cat("update", t, ":", s, "with", a, ": Q", signif(Q[s, a], 3))
      
      Q[s, a] <-
        Q[s, a] + alpha_val * (r + gamma * max(Q[sp, ]) - Q[s, a])
    
      if (verbose > 1) 
        cat(" ->", signif(Q[s, a], 3), "- alpha:", signif(alpha_val, 3), "\n")
      
    }
    
    # return is handled by on.exit()
  }
