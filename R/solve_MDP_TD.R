# Solve MDPs using Temporal Differencing

#' @rdname solve_MDP
#' @param alpha step size as a function of the time step `t` and the number of times
#'   the respective Q-value was updated `n` or a scalar. For expected Sarsa, alpha is 
#'   often set to 1.
#' @param epsilon used for \eqn{\epsilon}-greedy policies.
#' @param n number of episodes used for learning.
#' @param Q a state-action value matrix.
#' @export
solve_MDP_TD <-
  function(model,
           method = "q_learning",
           horizon = NULL,
           discount = NULL,
           alpha = function(t, n) min(10/n, 1),
           epsilon = 0.2,
           n = 1000,
           Q = 0,
           matrix = TRUE,
           continue = FALSE,
           progress = TRUE,
           verbose = FALSE,
           ...) {
    .nodots(...)
    
    method <-
      match.arg(method, c("sarsa", "q_learning", "expected_sarsa"))
    
    if (verbose > 1)
      progress <- FALSE
    
    if (is.null(horizon))
      horizon <- model$horizon
    if (is.null(horizon) || !is.finite(horizon)) {
      stop("Finite horizon needed for exploration!")
    }
    
    if (is.null(discount)) {
      discount <- model$discount
    }
    if (is.null(discount)) {
      discount <- 1
    }
    model$discount <- discount
   
    if (!is.function(alpha))
      alpha_val <- alpha
     
    if (matrix) {
      if (verbose)
        cat("Precomputing matrices for transitions and rewards ...")
      model <- normalize_MDP(
        model,
        sparse = NULL,
        precompute_absorbing = TRUE,
        progress = progress
      )
      
      if (verbose)
        cat(" done.\n")
    }
    
    if (progress) {
      pb <- my_progress_bar(n + 1L, name = "solve_MDP")
      pb$tick(0)
    }
      
    S <- S(model)
    A <- A(model)
    start <- start_vector(model, sparse = FALSE)

    # Initialize Q and Q_N
    if (continue) {
      if (is.null(model$solution$Q))
        stop("model solution does not contain a Q matrix to continue from!")
      Q <- model$solution$Q
      Q_N <- model$solution$Q_N
    } else {
      Q <- init_Q(model, Q)
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
      
      if (e < n)
          warning("Manual interupt: MDP solver stopped at episode ", e)
      
      if (verbose) {
        cat("\nTerminated at episode:", e, "\n")
      }
      
      model$solution <- list(
        method = method,
        alpha = alpha,
        epsilon = epsilon,
        n = n,
        Q = Q,
        Q_N = Q_N,
        converged = NA,
        policy = list(greedy_policy(Q))
      )
      return(model)
    })
    
    if (verbose) {
      cat("Running", method)
      cat("\nalpha:            ", deparse(alpha))
      cat("\nepsilon:          ", epsilon)
      cat("\nn                 ", n, "\n")
      
      cat("\nInitial Q (first 20 max):\n")
      print(head(Q, n = 20))
      
      cat("\nInitial Q_N (first 20 max):\n")
      print(head(Q_N, n = 20))
      cat("\n")
    }
    
    # loop episodes
    e <- 0L
    while (e < n) {
      e <- e +1L
      if (progress)
        pb$tick()
      
      s <- sample(S, 1L, prob = start)
      a <- greedy_action(Q, s, epsilon)

      # loop steps in episode
      i <- 0L
      while (TRUE) {
        i <- i + 1L
        if ((i %% 100 == 0) && progress && !pb$finished)
          pb$tick(0)
        
        #s_prime <- sample(S, 1L, prob = P[[a]][s, ])
        s_prime <- sample(S, 1L, prob = transition_matrix(model, a, s, sparse = FALSE))
        #s_prime <- sample_sparse(S, 1L, prob = transition_matrix(model, a, s, sparse = NULL))
        
        r <- reward_matrix(model, a, s, s_prime)
        a_prime <- greedy_action(Q, s_prime, epsilon)

        if (verbose > 1) {
          if (i == 1L) {
            cat("\n*** Episode", e, "***\n")
          }
          cat(sprintf("Ep %i step %i - transition: %s -> %s; action: %i -> %i q(s,a): %.3f -> ", 
                      e, i, s, s_prime, a, a_prime, Q[s, a]))
        }

        # update Q and Q_N and calculate alpha
        Q_N[s, a] <- Q_N[s, a] + 1L
        if (is.function(alpha))
          #alpha_val <- alpha(t, Q_N[s, a])
          alpha_val <- alpha(t, sum(Q_N[s, ]))
        
        if (method == "sarsa") {
          # on-policy: is called Sarsa because it uses the sequence s, a, r, s', a'
          Q[s, a] <-
            Q[s, a] + alpha_val * (r + discount * Q[s_prime, a_prime] - Q[s, a])
        
        } else if (method == "q_learning") {
          # we update the greedy action (max Q) -> off policy
          Q[s, a] <-
            Q[s, a] + alpha_val * (r + discount * max(Q[s_prime, ]) - Q[s, a])
        
        } else if (method == "expected_sarsa") {
          # on-policy 
          exp_Q_prime <-
            sum(greedy_action(Q, s_prime, epsilon, prob = TRUE) * Q[s_prime, ],
                na.rm = TRUE
            )
          Q[s, a] <-
            Q[s, a] + alpha_val * (r + discount * exp_Q_prime - Q[s, a])
        }

        if (is.na(Q[s, a])) {
          Q[s, a] <- -Inf
        }
        
        if (verbose > 1) {
          cat(sprintf("%.3f (N: %i alpha: %.3f)\n", Q[s, a], Q_N[s,a], alpha_val))
        }
        
        
        s <- s_prime
        a <- a_prime

        if (absorbing_states(model, state = s))
          break

        if (i >= horizon)
          break

      }
    }
    
    # return via on.exit()
  }
