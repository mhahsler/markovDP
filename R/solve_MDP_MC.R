# Solve MDPs using Monte Carlo Methods

if (FALSE) {
  data(Maze)
  sol_def <- solve_MDP(Maze)
  policy(sol_def)
  gridworld_plot(sol_def)
  
  sol <- solve_MDP_MC(Maze,
                      horizon = 100,
                      N = 1000,
                      first_visit = TRUE,
                      verbose = FALSE)
  policy(sol)
  gridworld_plot(sol)
  
  sum(abs(policy(sol_def)$U - policy(sol)$U))
  
  sol <- solve_MDP_MC(Maze,
                      method = "MC_on_policy",
                      horizon = 100,
                      N = 1000,
                      first_visit = TRUE,
                      verbose = FALSE)
  policy(sol)
  gridworld_plot(sol)
  sum(abs(policy(sol_def)$U - policy(sol)$U))
  
}

#' @rdname solve_MDP
#' @export
solve_MDP_MC <-
  function(model,
           method = "MC_exploring_starts",
           horizon = NULL,
           discount = NULL,
           N = 100,
           ...,
           U = NULL,
           first_visit = TRUE,
           verbose = FALSE) {
    ### default is infinite horizon, but we use 10000 to guarantee termination
    warn_horizon <- FALSE
    if (is.null(horizon)) {
      horizon <- model$horizon
    }
    if (is.null(horizon)) {
      if (!any(absorbing_states(model))) {
        stop("The model has no absorbing states. Specify the horizon.")
      }
      warn_horizon <- TRUE
      horizon <- 10000
    }
    
    if (is.null(discount)) {
      discount <- model$discount
    }
    if (is.null(discount)) {
      discount <- 1
    }
    model$discount <- discount
    
    
    method <-
      match.arg(method,
                c("MC_exploring_starts", "MC_on_policy", "MC_off_policy"))
    
    
    switch(
      method,
      MC_exploring_starts = MC_exploring_starts(model, method, horizon, discount, N, U, first_visit, verbose),
      MC_on_policy = MC_exploring_starts(model, method, horizon, discount, N, U, first_visit, verbose, ...),
      MC_off_policy = stop("Not implemented")
    )
  }

MC_exploring_starts <- function(model,
                                method,
                                horizon,
                                discount,
                                N,
                                U = NULL,
                                first_visit = TRUE,
                                verbose = FALSE) {
  ## Learns a greedy policy. In order to still keep exploring it uses the
  ## idea of exploring starts: All state-action pairs have a non-zero
  ## probability of being selected as the start of an episode.
  ## (RL book, Chapter 5)
  
  S <- model$states
  A <- model$actions
  gamma <- discount
  
  # Start with arbitrary policy
  pi <- random_policy(model)
  
  if (is.null(U)) {
    Q <- matrix(0,
                nrow = length(S),
                ncol = length(A),
                dimnames = list(S, A))
  } else {
    Q <- q_values(model, U = U)
  }
  
  # instead of returns we use a more efficient running average where Q_N is
  # the number of averaged values.
  Q_N <- matrix(0L,
                nrow = length(S),
                ncol = length(A),
                dimnames = list(S, A))
  
  if (verbose) {
    cat("Initial policy:\n")
    print(pi)
    
    cat("Initial Q:\n")
    print(Q)
  }
  
  # Loop through N episodes
  for (e in seq(N)) {
    
    # add faster without checks
    #model <- add_policy(model, policy = pi)
    model$solution <- list(
      method = "manual",
      policy = list(pi),
      converged = NA
    )
    
    ep <- simulate_MDP(
      model,
      n = 1,
      horizon = horizon,
      return_trajectories = TRUE,
      exploring_starts = TRUE
    )$trajectories
    
    if (verbose) {
      cat(paste("*** Episode", e, "***\n"))
      print(ep)
    }
    
    G <- 0
    for (i in rev(seq_len(nrow(ep)))) {
      r_t_plus_1 <- ep$r[i]
      s_t <- ep$s[i]
      a_t <- ep$a[i]
      
      G <- gamma * G + r_t_plus_1
      
      # Only update for first visit of a s/a combination
      if (!first_visit || i < 2L ||
          !any(s_t == ep$s[1:(i - 1L)] &
               a_t == ep$a[1:(i - 1L)])) {
        if (verbose)
          cat(paste0(
            "Update at step ",
            i,
            ":\n",
            "  - Q(",
            s_t,
            ", ",
            a_t,
            "): ",
            round(Q[s_t, a_t], 3)
          ))
        
        # running average instead of averaging Returns lists.
        Q[s_t, a_t] <- (Q[s_t, a_t] * Q_N[s_t, a_t] + G) / (Q_N[s_t, a_t] +
                                                              1)
        Q_N[s_t, a_t] <- Q_N[s_t, a_t] + 1L
        
        if (verbose)
          cat(paste0(" -> ", round(Q[s_t, a_t], 3), " (G = ", round(G, 3), ")\n"))
        
        if (verbose)
          cat(paste0("  - pi[", s_t, "]: ", pi[s_t, ]))
        
        pi$action[s_t] <- greedy_action(Q, s_t)
        
        if (verbose)
          cat(paste0(" -> ", pi[s_t, ], "\n"))
      }
    }
  }
  
  model$solution <- list(
    method = method,
    N = N,
    Q = Q,
    converged = NA,
    policy = list(greedy_policy(Q))
  )
  
  model
}

MC_on_policy <- function(model,
                                method,
                                horizon,
                                discount,
                                N,
                                U = NULL,
                                first_visit = TRUE,
                                verbose = FALSE,
                                epsilon) {
  ## Learns an epsilon-soft policy (also used as behavior)
  ## (RL book, Chapter 5)
  
  S <- model$states
  A <- model$actions
  gamma <- discount
  
  # Start with arbitrary policy, we make it soft by specifying epsilon 
  # in the simulation
  pi <- random_policy(model)
  
  if (is.null(U)) {
    Q <- matrix(0,
                nrow = length(S),
                ncol = length(A),
                dimnames = list(S, A))
  } else {
    Q <- q_values(model, U = U)
  }
  
  # instead of returns we use a more efficient running average where Q_N is
  # the number of averaged values.
  Q_N <- matrix(0L,
                nrow = length(S),
                ncol = length(A),
                dimnames = list(S, A))
  
  if (verbose) {
    cat("Initial policy:\n")
    print(pi)
    
    cat("Initial Q:\n")
    print(Q)
  }
  
  # Loop through N episodes
  for (e in seq(N)) {
    
    # add faster without checks
    #model <- add_policy(model, policy = pi)
    model$solution <- list(
      method = "manual",
      policy = list(pi),
      converged = NA
    )
    
    # use epsilon-soft policy!
    ep <- simulate_MDP(
      model,
      n = 1,
      horizon = horizon,
      epsilon = epsilon,
      return_trajectories = TRUE
    )$trajectories
    
    if (verbose) {
      cat(paste("*** Episode", e, "***\n"))
      print(ep)
    }
    
    G <- 0
    for (i in rev(seq_len(nrow(ep)))) {
      r_t_plus_1 <- ep$r[i]
      s_t <- ep$s[i]
      a_t <- ep$a[i]
      
      G <- gamma * G + r_t_plus_1
      
      # Only update for first visit of a s/a combination
      if (!first_visit || i < 2L ||
          !any(s_t == ep$s[1:(i - 1L)] &
               a_t == ep$a[1:(i - 1L)])) {
        if (verbose)
          cat(paste0(
            "Update at step ",
            i,
            ":\n",
            "  - Q(",
            s_t,
            ", ",
            a_t,
            "): ",
            round(Q[s_t, a_t], 3)
          ))
        
        # running average instead of averaging Returns lists.
        Q[s_t, a_t] <- (Q[s_t, a_t] * Q_N[s_t, a_t] + G) / (Q_N[s_t, a_t] +
                                                              1)
        Q_N[s_t, a_t] <- Q_N[s_t, a_t] + 1L
        
        if (verbose)
          cat(paste0(" -> ", round(Q[s_t, a_t], 3), " (G = ", round(G, 3), ")\n"))
        
        if (verbose)
          cat(paste0("  - pi[", s_t, "]: ", pi[s_t, ]))
        
        # the simulation takes care of the epsilon
        pi$action[s_t] <- greedy_action(Q, s_t)
        
        if (verbose)
          cat(paste0(" -> ", pi[s_t, ], "\n"))
      }
    }
  }
  
  model$solution <- list(
    method = method,
    N = N,
    Q = Q,
    converged = NA,
    policy = list(greedy_policy(Q))
  )
  
  model
}