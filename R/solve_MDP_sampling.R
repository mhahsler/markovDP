# Solve MDPs using Random-Sample one-step tabular Q-planning
#(Sutton & Barto, Chapter 8)

#' @rdname solve_MDP
#' @details
#' ## Planning by Sampling
#'
#' A simple, not very effective, planning method proposed by Sutton and Barto (2020) in Chapter 8.
#'
#' * **Random-sample one-step tabular Q-planning** randomly selects a
#' state/action pair and samples the resulting reward and next state from
#' the model. This
#' information is used to update a single Q-table value.
#' 
#' @export
solve_MDP_sampling <-
  function(model,
           method = "q_planning",
           horizon = NULL,
           discount = NULL,
           alpha = function(t, n) min(10/n, 1),
           n = 1000,
           Q = NULL,
           matrix = TRUE,
           continue = FALSE,
           progress = TRUE,
           verbose = FALSE,
           ...) {
    .nodots(...)
    
    if (verbose > 1)
      progress <- FALSE
    
    if (matrix) {
      if (verbose)
        cat("Precomputing matrices for R and T ...")
      model <- normalize_MDP(
        model,
        sparse = NULL,
        precompute_absorbing = FALSE,
        progress = progress
      )
      
      if (verbose)
        cat(" done.\n")
    }
    
    PROGRESS_INTERVAL <- 100
    
    if (!is.function(alpha))
      alpha_val <- alpha
    
    if (!is.null(horizon))
      warning("q_planning does not use horizon. The specified horizon is ignored.")
    
    discount <- discount %||% model$discount %||% 1
    model$discount <- discount
    gamma <- discount
    
    S <- S(model)
    A <- A(model)
    start <- start_vector(model, sparse = FALSE)
    
    method <-
      match.arg(method, c("q_planning"))
    
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
    
    # loop through tries
    t <- 0L
    while (t < n) {
      t <- t + 1L
      if (!(t %% PROGRESS_INTERVAL) && progress)
        pb$tick()
      
      # sample a state/action pair
      s <- sample(S, 1L)
      a <- sample(A, 1L)
      
      sp_r <- act.int(model, s, a)
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
