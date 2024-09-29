# Solve MDPs using Monte Carlo Methods

#' @rdname solve_MDP
#' @param first_visit if `TRUE` then only the first visit of a state/action pair
#'   in an episode is used to update Q, otherwise, every-visit update is used. 
#' @export
solve_MDP_MC <-
  function(model,
           method = "MC_exploring_starts",
           horizon = NULL,
           discount = NULL,
           n = 100,
           ...,
           Q = NULL,
           epsilon = 0.1,
           first_visit = TRUE,
           continue = FALSE,
           progress = TRUE,
           verbose = FALSE) {
    method <-
      match.arg(method,
                c("MC_exploring_starts", "MC_on_policy", "MC_off_policy"))
 
    ### default is infinite horizon, but we use 1000 to guarantee termination
    if (is.null(horizon)) {
      horizon <- model$horizon
    }
    if (is.null(horizon) || is.infinite(horizon)) {
      warning("No finite horizon defined. Using a maximum horizon of 1000 to guarantee termination. Specify horizon to remove this warning.")
      model$horizon <- horizon <- 1000
    }
    
    if (is.null(discount)) {
      discount <- model$discount
    }
    if (is.null(discount)) {
      discount <- 1
    }
    model$discount <- discount
    
    # Initialize Q
    if (continue) {
      if (is.null(model$solution$Q))
        stop("model solution does not contain a Q matrix to continue from!")
      Q <- model$solution$Q
    } else if (is.null(Q)) {
      Q <-
        matrix(0,
               nrow = length(model$states),
               ncol = length(model$actions),
               dimnames = list(model$states, model$actions)
        )
    } 
    
    switch(
      method,
      MC_exploring_starts = MC_on_policy(model, method, horizon, 
                                  discount, n, Q, exploring_starts = TRUE,
                                  epsilon = 0, first_visit = first_visit, 
                                  progress = progress, verbose = verbose, ...),
      MC_on_policy = MC_on_policy(model, method, horizon, 
                                  discount, n, Q, exploring_starts = FALSE,
                                  epsilon = epsilon, first_visit = first_visit, 
                                  progress = progress, verbose = verbose, ...),
      MC_off_policy = MC_off_policy(model, method, horizon, 
                                  discount, n, Q, first_visit, 
                                  progress, verbose, epsilon = epsilon, ...)
    )
  }

MC_on_policy <- function(model,
                                method,
                                horizon,
                                discount,
                                n,
                                Q = NULL,
                                exploring_starts,
                                epsilon,
                                first_visit = TRUE,
                                progress = TRUE,
                                verbose = FALSE) {
  ## exploring starts: Learns a greedy policy. In order to still keep exploring it uses the
  ## idea of exploring starts: All state-action pairs have a non-zero
  ## probability of being selected as the start of an episode.
  ## (RL book, Chapter 5)
  
  ## on policy: Learns an epsilon-soft policy (also used as behavior)
  ## (RL book, Chapter 5)
  
  if (exploring_starts) {
    if (epsilon != 0)
      warning("epsilon should be 0 for exploring starts!")
  }
  
  S <- model$states
  A <- model$actions
  gamma <- discount
  
  # Start with arbitrary policy, we make it soft by specifying epsilon 
  # in the simulation
  
  if (is.null(Q)) {
    Q <- matrix(0,
                nrow = length(S),
                ncol = length(A),
                dimnames = list(S, A))
    pi <- random_policy(model, only_available_actions = TRUE)
  } else {
    pi <- greedy_policy(Q)
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
  
  if (progress)
    pb <- my_progress_bar(n, name = "solve_MDP")
  
  on.exit({ 
    warning("MDP solver manually interrupted early.")
    
    if (verbose) {
      cat("\nTerminated during episode:", e, "\n")
    }
    
    model$solution <- list(
      method = method,
      n = n,
      Q = Q,
      converged = NA,
      policy = list(greedy_policy(Q))
    )
    
    return(model)
  })
  
  # Loop through N episodes
  for (e in seq(n)) {
    if (progress)
      pb$tick()
    
    # add faster without checks
    #model <- add_policy(model, policy = pi)
    model$solution <- list(
      method = "manual",
      policy = list(pi),
      converged = NA
    )
    
    # use epsilon-soft policy!
    ep <- sample_MDP(
      model,
      n = 1,
      horizon = horizon,
      epsilon = epsilon,
      exploring_starts = exploring_starts,
      return_trajectories = TRUE,
      progress = FALSE,
      verbose = FALSE
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
          cat(paste0("  - pi(", s_t, "): ", pi[s_t, "action"]))
        
        # the simulation takes care of the epsilon
        pi$action[s_t] <- greedy_action(Q, s_t)
        
        if (verbose)
          cat(paste0(" -> ", pi[s_t, "action"], "\n"))
      }
    }
  }
  
  on.exit()
  
  model$solution <- list(
    method = method,
    n = n,
    Q = Q,
    converged = NA,
    policy = list(greedy_policy(Q))
  )
  
  model
}


MC_off_policy <- function(model,
                          method,
                          horizon,
                          discount,
                          n,
                          Q = NULL,
                          first_visit = TRUE,
                          progress = TRUE,
                          verbose = FALSE,
                          epsilon) {
  ## Learns an epsilon-greedy policy using an epsilon-soft policy for behavior
  ## (RL book, Chapter 5)
  
  S <- model$states
  A <- model$actions
  gamma <- discount
  
  # Initialize
  if (is.null(Q)) {
    Q <- matrix(runif(length(S) * length(A)), nrow = length(S), ncol = length(A), dimnames = list(S, A))
  }
  
  pi <- greedy_policy(Q)
  
  # cumulative sum of the weights W used in incremental updates
  C <- matrix(0L, nrow = length(S), ncol = length(A), dimnames = list(S, A))
  
  if (progress)
    pb <- my_progress_bar(n, name = "solve_MDP")
  
  on.exit({ 
    warning("MDP solver manually interrupted early.")
    
    if (verbose) {
      cat("\nTerminated during episode:", e, "\n")
    }
    
    model$solution <- list(
      method = method,
      n = n,
      Q = Q,
      converged = NA,
      policy = list(greedy_policy(Q))
    )
    
    return(model)
  })
  
  # Loop through N episodes
  for (e in seq(n)) {
    if (progress)
      pb$tick()
    
    # we use as the soft behavioral policy an epsilon-soft version of pi.
    # use epsilon-soft policy!
    b <- pi
    
    # add faster without checks
    #model <- add_policy(model, policy = b)
    model$solution <- list(
      method = "manual",
      policy = list(b),
      converged = NA
    )
    
    ep <- sample_MDP(
      model,
      n = 1,
      horizon = horizon,
      epsilon = epsilon,
      return_trajectories = TRUE,
      progress = FALSE,
      verbose = FALSE
    )$trajectories
     
    if (verbose) {
      cat(paste("*** Episode", e, "***\n"))
      print(ep)
    }
    
    G <- 0
    W <- 1
    for (i in rev(seq_len(nrow(ep)))) {
      r_t_plus_1 <- ep$r[i]
      s_t <- ep$s[i]
      a_t <- ep$a[i]
      
      G <- gamma * G + r_t_plus_1
      
      # increase cumulative sum of W and update Q with weighted G
      C[s_t, a_t] <- C[s_t, a_t] + W
      Q[s_t, a_t] <- Q[s_t, a_t] + (W/C[s_t, a_t]) * (G - Q[s_t, a_t])
      
      pi$action[s_t] <- greedy_action(Q, s_t)
      
      # the algorithm can only learn from the tail of the episode where b
      # also used the greedy actions in pi. The method is inefficient and
      # cannot use all the data!
      if (a_t != pi$action[s_t])
        break
      
      # update the weight using b(A_t|S_t)
      # Note, we could used available_actions(model, s_t), but that is expensive 
      if (a_t == b$action[s_t]) 
        b_at_st <- 1 - epsilon + epsilon / length(A)
      else
        b_at_st <- epsilon / length(A)
       
      
      W <- W * 1/b_at_st
    }
  }
  
  on.exit()
  
  model$solution <- list(
    method = method,
    n = n,
    Q = Q,
    converged = NA,
    policy = list(greedy_policy(Q))
  )
  
  model
  }