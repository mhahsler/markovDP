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
           epsilon = NULL,
           alpha = NULL,
           first_visit = TRUE,
           matrix = TRUE,
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
    } 
    
    if (matrix) {
      if (verbose)
        cat("Precomputing dense matrices for R and T ...")
      model <- normalize_MDP(
        model,
        sparse = FALSE,
        precompute_absorbing = TRUE,
        progress = progress
      )
      
      if (verbose)
        cat(" done.\n")
    }
    
    switch(
      method,
      MC_exploring_starts = MC_on_policy(model, method, horizon, 
                                  discount, n, Q, exploring_starts = TRUE,
                                  epsilon = epsilon, alpha = alpha, first_visit = first_visit, 
                                  progress = progress, verbose = verbose, ...),
      MC_on_policy = MC_on_policy(model, method, horizon, 
                                  discount, n, Q, exploring_starts = FALSE,
                                  epsilon = epsilon, alpha = alpha, first_visit = first_visit, 
                                  progress = progress, verbose = verbose, ...),
      MC_off_policy = MC_off_policy(model, method, horizon, 
                                  discount, n, Q, 
                                  epsilon = epsilon, alpha = alpha, first_visit = first_visit, 
                                  progress = progress, verbose = verbose, ...)
    )
  }

MC_on_policy <- function(model,
                                method,
                                horizon,
                                discount,
                                n,
                                Q = NULL,
                                exploring_starts,
                                epsilon = NULL,
                                alpha = NULL,
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
    epsilon <- epsilon %||% 0
    if (epsilon != 0)
      warning("epsilon should be 0 for exploring starts!")
  } else {
    epsilon <- epsilon %||% .2
  }
  
  S <- S(model)
  A <- A(model)
 
  alpha_param <- alpha
   
  # Start with arbitrary policy, we make it soft by specifying epsilon 
  # in the simulation.
  # Instead of returns we use a more efficient running average where Q_N is
  # the number of averaged values.
  
  if (is.null(Q)) {
    Q <- Q_zero(model)
    Q_N <-
      matrix(0L,
             nrow = nrow(Q),
             ncol = ncol(Q),
             dimnames = dimnames(Q)
      )
    pi <- random_policy(model, only_available_actions = TRUE)
  } else {
    pi <- greedy_policy(Q)
    Q_N <- model$solution$Q_N
    if (is.null(Q_N))
      stop("Q_N missing in previous solution. Cannot continue!")
  }
  
  if (verbose) {
    cat("Running MC_on_policy\n")
    cat("\nalpha:            ", alpha %||% "1/N")
    cat("\nepsilon:          ", epsilon)
    cat("\nexploring starts: ", exploring_starts, "\n")
    
    cat("\nInitial policy (first 20 max):\n")
    print(head(pi, n = 20))
    
    cat("\nInitial Q (first 20 max):\n")
    print(head(Q, n = 20))
    
    cat("\nInitial Q_N (first 20 max):\n")
    print(head(Q_N, n = 20))
    cat("\n")
  }
  
  if (progress)
    pb <- my_progress_bar(n + 1L, name = "solve_MDP")
  
  on.exit({
    if (progress) {
      pb$tick(0)
      pb$terminate()
    }
    
    if (e < n)
      warning("Manual interupt: MDP solver stopped at episode ", e)
    
    if (verbose) {
      cat("\nTerminated after episode:", e, "\n")
    }
    
    model$solution <- list(
      method = method,
      n = n,
      Q = Q,
      Q_N = Q_N,
      converged = NA,
      policy = list(greedy_policy(Q))
    )
    
    return(model)
  })
  
  # Loop through N episodes
  e <- 0L
  while (e < n) {
    e <- e + 1L
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
      trajectories = TRUE,
      progress = FALSE,
      verbose = FALSE
    )$trajectories
    
    if (verbose) {
      cat(paste("\n****************** Episode", e, "******************\n"))
      print(ep)
      cat("\n")
    }
    
    G <- 0
    for (i in rev(seq_len(nrow(ep)))) {
      r_t_plus_1 <- ep$r[i]
      s_t <- ep$s[i]
      a_t <- ep$a[i]
      
      G <- discount * G + r_t_plus_1
      
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
        # Q <- Q + 1/n (G - Q) ... 1/n is alpha
        Q_N[s_t, a_t] <- Q_N[s_t, a_t] + 1L
        
        if (is.null(alpha_param))
          alpha <- 1 / Q_N[s_t, a_t]
            
        err <- G - Q[s_t, a_t]
        if (!is.nan(err))
          Q[s_t, a_t] <- Q[s_t, a_t] + alpha * (err)
        
        if (verbose)
          cat(paste0(" -> ", round(Q[s_t, a_t], 3), 
                     " (G = ", round(G, 3), 
                     "; alpha = ", signif(alpha,3), ")\n"))
        
        if (verbose)
          cat(paste0("  - pi(", s_t, "): ", pi[s_t, "action"]))
        
        # the simulation takes care of the epsilon
        pi$action[s_t] <- greedy_action(Q, s_t)
        
        if (verbose)
          cat(paste0(" -> ", pi[s_t, "action"], "\n"))
      }
    }
  }
  
  # return is handled by on.exit()
}


MC_off_policy <- function(model,
                          method,
                          horizon,
                          discount,
                          n,
                          Q = NULL,
                          epsilon = NULL,
                          alpha = NULL,
                          first_visit = TRUE,
                          progress = TRUE,
                          verbose = FALSE,
                          ...
                          ) {
  ## Learns an epsilon-greedy policy using an epsilon-soft policy for behavior
  ## (RL book, Chapter 5)
  
  if (!is.null(alpha))
    warning("MC_off_policy does not use alpha.")
  
  epsilon <- epsilon %||% .2
  
  S <- S(model)
  A <- A(model)
  gamma <- discount
  
  # Initialize
  if (is.null(Q)) {
    Q <- Q_zero(model)
    C <- matrix(0L, nrow = nrow(Q), ncol = ncol(Q), dimnames = dimnames(Q))
  } else { # we get Q for continuation
    C <- model$solution$C
    if (is.null(C))
      stop("C missing in previous solution. Cannot continue!")
  }
  
  pi <- greedy_policy(Q)
  
  # cumulative sum of the weights W used in incremental updates
  
  if (verbose) {
    cat("Running MC_off_policy\n")
    
    cat("\nepsilon: ", epsilon, "\n")
    
    cat("\nInitial policy (first 20 max):\n")
    print(head(pi, n = 20))
    
    cat("\nInitial Q (first 20 max):\n")
    print(head(Q, n = 20))
    
    cat("\nInitial C (first 20 max):\n")
    print(head(C, n = 20))
    cat("\n")
  }
  
  
  if (progress)
    pb <- my_progress_bar(n + 1L, name = "solve_MDP")
  
  on.exit({ 
    if (progress) {
      pb$tick(0)
      pb$terminate()
    }
    
    if (e < n)
      warning("Manual interupt: MDP solver stopped at episode ", e)
    
    if (verbose) {
      cat("\nTerminated after episode:", e, "\n")
    }
    
    model$solution <- list(
      method = method,
      n = n,
      Q = Q,
      C = C,
      converged = NA,
      policy = list(greedy_policy(Q))
    )
    
    return(model)
  })
  
  # Loop through episodes
  e <- 0L
  while (e < n) {
    e <- e + 1L
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
      trajectories = TRUE,
      progress = FALSE,
      verbose = FALSE
    )$trajectories
     
    if (verbose) {
      cat(paste("\n****************** Episode", e, "******************\n"))
      print(ep)
      cat("\n")
    }
    
    G <- 0
    W <- 1
    
    for (i in rev(seq_len(nrow(ep)))) {
      
      r_t_plus_1 <- ep$r[i]
      s_t <- ep$s[i]
      a_t <- ep$a[i]
        
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
      
      G <- gamma * G + r_t_plus_1
      
      # increase cumulative sum of W and update Q with weighted G
      C[s_t, a_t] <- C[s_t, a_t] + W
      Q[s_t, a_t] <- Q[s_t, a_t] + (W / C[s_t, a_t]) * (G - Q[s_t, a_t])
        
      if (verbose)
        cat(paste0(" -> ", round(Q[s_t, a_t], 3), 
                   " (G = ", round(G, 3), 
                   "; W/C = ", signif(W / C[s_t, a_t], 3), ")\n"))
      
      if (verbose)
        cat(paste0("  - pi(", s_t, "): ", pi[s_t, "action"]))
      
      pi$action[s_t] <- greedy_action(Q, s_t)
      
      if (verbose)
        cat(paste0(" -> ", pi[s_t, "action"], "\n"))
      
      # the algorithm can only learn from the tail of the episode where b
      # also used the greedy actions in pi. The method is inefficient and
      # cannot use all the data!
      if (a_t != pi$action[s_t]) {
        if (verbose)
          cat("Break: a_t is not the greedy action.\n")
        break
      }
      
      # update the weight W = pi(A_t|S_t)/b(A_t|S_t) using pi = 1 and  b(A_t|S_t)
      # Note, we could used available_actions(model, s_t), but that is expensive 
      if (a_t == b$action[s_t]) 
        b_at_st <- 1 - epsilon + epsilon / length(A)
      else
        b_at_st <- epsilon / length(A)
       
      
      W <- W * 1 / b_at_st
    }
  }
  
  # return is handled by on.exit()
  }