#' Solve MDPs using Monte Carlo Control
#'
#' Solve MDPs using Monte Carlo control.
#'
#' @details
#' The idea is to estimate the action value function for a policy as the 
#' average of sampled returns.
#' 
#' \deqn{q_\pi(s,a) = \mathbb{E}_\pi[R_i|S_0=s,A_0=a] \approx \frac{1}{n} \sum_{i=1}^n R_i}
#' 
#' Monte Carlo control simulates a whole episode using the current behavior
#' policy and uses the sampled reward to update the Q values. For on-policy 
#' methods, the behavior policy is updated to be greedy (i.e., optimal) with 
#' respect to the new Q values. Then the next episode is simulated till 
#' the predefined number of episodes is completed.  
#' 
#' ## Implemented methods
#' 
#' Implemented are the following temporal difference control methods
#' described in Sutton and Barto (2018).
#'
#' * **Monte Carlo Control with exploring Starts** learns the optimal greedy policy. 
#' It uses the same greedy policy for
#' behavior and target (on-policy learning).
#' After each episode, the policy is updated to be greedy with respect to the 
#' current Q values. 
#' To make sure all states/action pairs are
#' explored, it uses exploring starts meaning that new episodes are started at a randomly
#' chosen state using a randomly chooses action.
#'
#' * **On-policy Monte Carlo Control** learns an epsilon-greedy policy
#' which it uses for behavior and as the target policy
#' (on-policy learning). An epsilon-greedy policy is used to provide 
#' exploration. For calculating running averages, an update with \eqn{\alpha = 1/n}
#' is used by default. A different update factor can be set using the parameter `alpha`
#' as either a fixed value or a function with the signature `function(t, n)` 
#' which returns the factor in the range \eqn{[0,1]}.
#'
#' * **Off-policy Monte Carlo Control** uses for behavior an arbitrary soft policy
#' (a soft policy has in each state a probability greater than 0 for all 
#' possible actions). 
#' We use an epsilon-greedy policy and the method learns a greedy policy using
#' importance sampling. Note: This method can only learn from the tail of the 
#' sampled runs where greedy actions are chosen. This means that it is very
#' inefficient in learning the beginning portion of long episodes. This problem 
#' is especially problematic when larger values for \eqn{\epsilon} are used. 
#'
#' #' ## Schedules
#' 
#' * epsilon schedule: `t` is increased by each processed episode.
#' * alpha schedule: `t` is set to the number of times the a Q-value for state-action
#'      combination was updated. 
#' 
#' @references 
#' 
#' Sutton, Richard S., and Andrew G. Barto. 2018. Reinforcement Learning: An Introduction. Second. The MIT Press. [http://incompleteideas.net/book/the-book-2nd.html](http://incompleteideas.net/book/the-book-2nd.html).
#' 
#' @inheritParams solve_MDP_TD
#' @param method string; one of the following solution methods:
#'    * `'exploring_starts'` - on-policy MC control with exploring starts.
#'    * `'on_policy'` - on-policy MC control with an \eqn{\epsilon}-greedy policy.
#'    * `'off_policy"'` - off-policy MC control using an \eqn{\epsilon}-greedy behavior policy.
#' @param first_visit if `TRUE` then only the first visit of a state/action pair
#'   in an episode is used to update Q, otherwise, every-visit update is used.
#' 
#' @inherit solve_MDP return 
#' 
#' @export
solve_MDP_MC <-
  function(model,
           method = "exploring_starts",
           horizon = NULL,
           discount = NULL,
           n = 100,
           Q = NULL,
           epsilon = NULL,
           alpha = NULL,
           first_visit = TRUE,
           ...,
           matrix = TRUE,
           continue = FALSE,
           progress = TRUE,
           verbose = FALSE) {
    .nodots(...)
    
    # this works for MDP and MDPTF with a defined state space
    if (!inherits(model, "MDPE") || is.null(S(model)))
      stop("The model needs to be an MDP description with a specified state space.")
    
    methods <- c("exploring_starts", "on_policy", "off_policy")
    method <- match.arg(method, methods)
    
    model <- .prep_model(model, horizon, discount, matrix, verbose, progress)
    
    if (is.infinite(model$horizon)) {
      warning(
        "No finite horizon defined. Specify horizon to remove this warning."
      )
      model$horizon <- convergence_horizon(model, n_updates = 1)
    }
    
    # Initialize Q
    if (continue) {
      if (is.null(model$solution$Q))
        stop("model solution does not contain a Q matrix to continue from!")
      Q <- model$solution$Q
    }
    
    switch(
      method,
      exploring_starts = solve_MDP_MC_on_policy(
        model,
        method,
        n,
        Q,
        exploring_starts = TRUE,
        epsilon = epsilon,
        alpha = alpha,
        first_visit = first_visit,
        progress = progress,
        verbose = verbose,
        ...
      ),
      on_policy = solve_MDP_MC_on_policy(
        model,
        method,
        n,
        Q,
        exploring_starts = FALSE,
        epsilon = epsilon,
        alpha = alpha,
        first_visit = first_visit,
        progress = progress,
        verbose = verbose,
        ...
      ),
      off_policy = solve_MDP_MC_off_policy(
        model,
        method,
        n,
        Q,
        epsilon = epsilon,
        alpha = alpha,
        first_visit = first_visit,
        progress = progress,
        verbose = verbose,
        ...
      )
    )
  }

solve_MDP_MC_on_policy <- function(model,
                         method,
                         n,
                         Q = NULL,
                         exploring_starts,
                         epsilon = NULL,
                         alpha = NULL,
                         first_visit = TRUE,
                         progress = TRUE,
                         verbose = FALSE,
                         ...) {
  .nodots(...)
  
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
    epsilon <- epsilon %||% schedule_exp(1, -log(1e-5)/n)
  }
  
  
  alpha <- alpha %||% schedule_exp(.2, .001)
  
  ### alpha/epsilon func
  alpha_func <- NULL
  if (is.function(alpha))
    alpha_func <- alpha
  
  epsilon_seq <- NULL
  if (is.function(epsilon))
    epsilon_seq <- epsilon(seq(1, n))
  else if(length(epsilon) == n)
    epsilon_seq <- epsilon
  ###
  
  S <- S(model)
  A <- A(model)
  horizon <- model$horizon
  discount <- model$discount
  
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
             dimnames = dimnames(Q))
    pi <- random_policy(model, only_available_actions = TRUE)
  } else {
    pi <- greedy_policy(Q)
    Q_N <- model$solution$Q_N
    if (is.null(Q_N))
      stop("Q_N missing in previous solution. Cannot continue!")
  }
  
  if (verbose) {
    cat("Running MC_on_policy")
    cat("\nalpha:            ", show_schedule(alpha)) 
    cat("\nepsilon:          ", show_schedule(epsilon))
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
    
    ### alpha/epsilon func
    if (!is.null(epsilon_seq)) 
      epsilon <- epsilon_seq[e]
    
    # alpha uses the count and not the episode number!
    #if (!is.null(alpha_seq)) 
    #  alpha <- alpha_seq[e]
    ###
    
    # add faster without checks
    #model <- add_policy(model, policy = pi)
    model$solution <- list(method = "manual",
                           policy = list(pi),
                           converged = NA)
    
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
    
    if (verbose > 1) {
      cat(paste(
        "\n****************** Episode",
        e,
        "******************\n"
      ))
      print(ep)
      cat("\n")
    }
    
    G <- 0
    for (i in rev(seq_len(nrow(ep)))) {
      r_t_plus_1 <- ep$r[i]
      s_t <- normalize_state_id(ep$s[i], model)
      a_t <- ep$a[i]
      
      G <- discount * G + r_t_plus_1
      
      # Only update for first visit of a s/a combination
      if (first_visit && i < 2L &&
          (any(s_t == ep$s[1:(i - 1L)] &
               a_t == ep$a[1:(i - 1L)])))
        next
      
      if (verbose > 1)
        cat(paste0(
          "Update at step ",
          i,
          ":",
          " Q(",
          s_t,
          ", ",
          a_t,
          "): ",
          round(Q[s_t, a_t], 3)
        ))
      
      # Q <- avg(Returns)
      # running average instead of averaging Returns lists.
      # Q <- Q + alpha (G - Q) ... default alpha is 1/n
      Q_N[s_t, a_t] <- Q_N[s_t, a_t] + 1L
      
      ### alpha/epsilon func
      # alpha uses the count and not the episode number!
      if (!is.null(alpha_func)) 
        alpha <- alpha_func(Q_N[s_t, a_t])
      ###
      
      err <- G - Q[s_t, a_t]
      if (!is.nan(err))
        Q[s_t, a_t] <- Q[s_t, a_t] + alpha * (err)
      
      if (verbose > 1)
        cat(paste0(
          " -> ",
          round(Q[s_t, a_t], 3),
          " (G = ",
          round(G, 3),
          "; alpha = ",
          signif(alpha, 3),
          ")"
        ))
      
      if (verbose > 1)
        cat(paste0("; pi(", s_t, "): ", pi[s_t, "action"]))
      
      # the simulation takes care of the epsilon
      pi$action[s_t] <- greedy_action(model, s_t, Q)
      
      if (verbose > 1)
        cat(paste0(" -> ", pi[s_t, "action"], "\n"))
      
    }
  }
  
  # return is handled by on.exit()
}


# Note: this is different from off-policy learning for Sarsa and does not use 
# the importance sampling ratio. It only uses the from the end of the episode as
# long as all actions match the greedy actions.

solve_MDP_MC_off_policy <- function(model,
                          method,
                          n,
                          Q = NULL,
                          epsilon = NULL,
                          alpha = NULL,
                          first_visit = TRUE,
                          progress = TRUE,
                          verbose = FALSE,
                          ...) {
  .nodots(...)
  
  # Learns an epsilon-greedy policy using an epsilon-soft policy for behavior
  ## (RL book, Chapter 5)
  
  if (!is.null(alpha))
    warning("MC_off_policy does not use alpha.")
  
  epsilon <- epsilon %||% .2
  
  S <- S(model)
  A <- A(model)
  gamma <- model$discount
  horizon <- model$horizon
  
  # Initialize
  if (is.null(Q)) {
    Q <- Q_zero(model)
    C <- matrix(0L,
                nrow = nrow(Q),
                ncol = ncol(Q),
                dimnames = dimnames(Q))
  } else {
    # we get Q for continuation
    C <- model$solution$C
    if (is.null(C))
      stop("C missing in previous solution. Cannot continue!")
  }
  
  pi <- greedy_policy(Q)
  
  # cumulative sum of the weights W used in incremental updates
  
  if (verbose) {
    cat("Running MC_off_policy")
    
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
    model$solution <- list(method = "manual",
                           policy = list(b),
                           converged = NA)
    
    ep <- sample_MDP(
      model,
      n = 1,
      horizon = horizon,
      epsilon = epsilon,
      trajectories = TRUE,
      progress = FALSE,
      verbose = FALSE
    )$trajectories
    
    if (verbose > 1) {
      cat(paste(
        "\n****************** Episode",
        e,
        "******************\n"
      ))
      print(ep)
      cat("\n")
    }
    
    G <- 0
    W <- 1
    
    for (i in rev(seq_len(nrow(ep)))) {
      r_t_plus_1 <- ep$r[i]
      s_t <- ep$s[i]
      a_t <- ep$a[i]
      
      if (verbose > 1)
        cat(paste0(
          "Update at step ",
          i,
          ":",
          " Q(",
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
      
      if (verbose > 1)
        cat(paste0(
          " -> ",
          round(Q[s_t, a_t], 3),
          " (G = ",
          round(G, 3),
          "; W/C = ",
          signif(W / C[s_t, a_t], 3),
          ");"
        ))
      
      if (verbose > 1)
        cat(paste0(" pi(", s_t, "): ", pi[s_t, "action"]))
      
      pi$action[s_t] <- greedy_action(model, s_t, Q)
      
      if (verbose > 1)
        cat(paste0(" -> ", pi[s_t, "action"], "\n"))
      
      # the algorithm can only learn from the tail of the episode where b
      # also used the greedy actions in pi. The method is inefficient and
      # cannot use all the data!
      if (a_t != pi$action[s_t]) {
        if (verbose > 1)
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