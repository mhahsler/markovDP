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
           N = 10000,
           U = NULL,
           progress = TRUE,
           verbose = FALSE
           ) {
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
    S_absorbing <- S[which(absorbing_states(model))]
    A <- model$actions
    start <- start_vector(model, sparse = FALSE)

    method <-
      match.arg(method, c("q_planning"))

    # Initialize Q
    if (is.null(U)) {
      Q <-
        matrix(0,
          nrow = length(S),
          ncol = length(A),
          dimnames = list(S, A)
        )
    } else {
      Q <- q_values(model, U = U)
    }

    if (progress)
      pb <- my_progress_bar(N, name = "solve_MDP")
    
    # loop through tries
    n <- N
    while (n > 0) {
      if (progress)
        pb$tick()
      
      n <- n - 1L
       
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

    model$solution <- list(
      method = method,
      alpha = alpha,
      N = N,
      Q = Q,
      converged = NA,
      policy = list(data.frame(
        state = S,
        U = apply(Q, MARGIN = 1, max),
        action = A[apply(Q, MARGIN = 1, which.max.random)],
        row.names = NULL
      ))
    )

    model
  }
