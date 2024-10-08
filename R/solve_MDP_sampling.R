# Solve MDPs using Random-Sample one-step tabular Q-planning
#(Sutton & Barto, Chapter 8)

#' @rdname solve_MDP
#' @export
solve_MDP_sampling <-
  function(model,
           method = "q_planning",
           horizon = NULL,
           discount = NULL,
           alpha = function(t, n) 1/n,
           n = 1000,
           Q = NULL,
           continue = FALSE,
           progress = TRUE,
           verbose = FALSE) {
    if (verbose)
      progress <- FALSE
    
    PROGRESS_INTERVAL <- 100
    
    if (!is.function(alpha))
      alpha_val <- alpha
    
    if (!is.null(horizon))
      warning("q_planning does not use horizon. The specified horizon is ignored.")
    
    if (is.null(discount)) {
      discount <- model$discount
    }
    if (is.null(discount)) {
      discount <- 1
    }
    gamma <- discount
    model$discount <- discount
    
    S <- model$states
    A <- model$actions
    start <- start_vector(model, sparse = FALSE)
    
    method <-
      match.arg(method, c("q_planning"))
    
    
    # Initialize Q
    if (continue) {
      if (is.null(model$solution$Q) || is.null(model$solution$N))
        stop("model solution does not contain a Q matrix or the N count matrix to continue from!")
      Q <- model$solution$Q
      N <- model$solution$N
    } else if (is.null(Q)) {
      Q <-
        matrix(0,
               nrow = length(S),
               ncol = length(A),
               dimnames = list(S, A)
        )
      N <-
        matrix(0L,
               nrow = length(S),
               ncol = length(A),
               dimnames = list(S, A)
        )
    }
    
    # return unconverged result when interrupted
    on.exit({
      warning("MDP solver interrupted early.")
      
      if (verbose) {
        cat("\nTerminated during iteration:", t, "\n")
      }
      
      model$solution <- list(
        method = method,
        alpha = alpha,
        n = n,
        Q = Q,
        N = N,
        converged = NA,
        policy = list(greedy_policy(Q))
      )
      return(model)
    })
    
    if (progress)
      pb <- my_progress_bar(n/PROGRESS_INTERVAL, name = "solve_MDP")
    
    # loop through tries
    for (t in seq(n)) {
      if (!(t %% PROGRESS_INTERVAL) && progress)
        pb$tick()
      
      # sample a state/action pair
      s <- sample(S, 1L)
      a <- sample(A, 1L)
      
      sp_r <- act(model, s, a)
      sp <- sp_r$state
      r <- sp_r$r
      
      # NOTE: we use as the default alpha = 1/N(s,a)
      #       then the expected error for each Q value is sigma_TD/sqrt(N(s,a))
      #           sigma_TD is the sd of R+gamma max_a Q(s', a)
     
      
      # update Q and N
      N[s, a] <- N[s, a] + 1L
      
      if (is.function(alpha))
        alpha_val <- alpha(t, N[s, a])
      
      if (verbose) 
        cat(s, "with", a, "(alpha:", signif(alpha_val, 3), ") Q:", signif(Q[s, a], 3))
      
      Q[s, a] <-
        Q[s, a] + alpha_val * (r + gamma * max(Q[sp, ]) - Q[s, a])
    
      if (verbose) 
        cat(" ->", signif(Q[s, a], 3), "\n")
      
    }
    
    on.exit()
    
    model$solution <- list(
      method = method,
      alpha = alpha,
      n = n,
      Q = Q,
      N = N,
      converged = NA,
      policy = list(greedy_policy(Q))
    )
    
    model
  }
