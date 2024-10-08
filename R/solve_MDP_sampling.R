# Solve MDPs using Random-Sample one-step tabular Q-planning
#(Sutton & Barto, Chapter 8)

#' @rdname solve_MDP
#' @export
solve_MDP_sampling <-
  function(model,
           method = "q_planning",
           horizon = NULL,
           discount = NULL,
           alpha = 0.5,
           n = 1000,
           Q = NULL,
           continue = FALSE,
           progress = TRUE,
           verbose = FALSE) {
    PROGRESS_INTERVAL <- 100
    
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
      if (is.null(model$solution$Q))
        stop("model solution does not contain a Q matrix to continue from!")
      Q <- model$solution$Q
    } else if (is.null(Q)) {
      Q <-
        matrix(0,
               nrow = length(S),
               ncol = length(A),
               dimnames = list(S, A)
        )
    }
    
    # return unconverged result when interrupted
    on.exit({
      warning("MDP solver interrupted early.")
      
      if (verbose) {
        cat("\nTerminated during iteration:", i, "\n")
      }
      
      model$solution <- list(
        method = method,
        alpha = alpha,
        n = n,
        Q = Q,
        converged = NA,
        policy = list(greedy_policy(Q))
      )
      return(model)
    })
    
    if (progress)
      pb <- my_progress_bar(n/PROGRESS_INTERVAL, name = "solve_MDP")
    
    # loop through tries
    for (i in seq(n)) {
      if (!(i %% PROGRESS_INTERVAL) && progress)
        pb$tick()
      
      # sample a state/act pair
      s <- sample(S, 1L)
      a <- sample(A, 1L)
      sp_r <- act(model, s, a)
      sp <- sp_r$state
      r <- sp_r$r
      
      # update Q
      Q[s, a] <-
        Q[s, a] + alpha * (r + gamma * max(Q[sp, , drop = FALSE]) - Q[s, a])
    }
    
    on.exit()
    
    model$solution <- list(
      method = method,
      alpha = alpha,
      n = n,
      Q = Q,
      converged = NA,
      policy = list(greedy_policy(Q))
    )
    
    model
  }
