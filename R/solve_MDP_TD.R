# Solve MDPs using Temporal Differencing

#' @rdname solve_MDP
#' @param alpha step size in `(0, 1]`.
#' @param epsilon used for \eqn{\epsilon}-greedy policies.
#' @param n number of episodes used for learning.
#' @param Q a action value matrix.
#' @export
solve_MDP_TD <-
  function(model,
           method = "q_learning",
           horizon = NULL,
           discount = NULL,
           alpha = 0.5,
           epsilon = 0.1,
           n = 1000,
           Q = NULL,
           continue = FALSE,
           progress = TRUE,
           verbose = FALSE) {
    method <-
      match.arg(method, c("sarsa", "q_learning", "expected_sarsa"))
    
    if (is.null(horizon))
      horizon <- model$horizon
    
    if (is.null(horizon) || !is.finite(horizon)) {
      horizon <- 100
      warning("Finite horizon needed, using ", horizon, immediate. = TRUE)
    }
    
    if (is.null(discount)) {
      discount <- model$discount
    }
    if (is.null(discount)) {
      discount <- 1
    }
    gamma <- discount
    model$discount <- discount
    
    if (progress) {
      pb <- my_progress_bar(n, name = "solve_MDP")
      pb$tick(0)
    }
      
    S <- model$states
    A <- model$actions
    #S_absorbing <- S[which(absorbing_states(model))]
    start <- start_vector(model, sparse = FALSE)


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
      warning("MDP solver manually interrupted early.")
      
      if (verbose) {
        cat("\nTerminated during iteration:", i, "\n")
      }
      
      model$solution <- list(
        method = method,
        alpha = alpha,
        epsilon = epsilon,
        n = n,
        Q = Q,
        converged = NA,
        policy = list(greedy_policy(Q))
      )
      return(model)
    })
    
    
    # loop episodes
    for (e in seq_len(n)) {
      if (progress)
        pb$tick()
      
      s <- sample(S, 1, prob = start)
      a <- greedy_action(Q, s, epsilon)

      # loop steps in episode
      i <- 1L
      while (TRUE) {
        if ((i %% 100 == 0) && !pb$finished && progress)
          pb$tick(0)
        #s_prime <- sample(S, 1L, prob = P[[a]][s, ])
        s_prime <- sample(S, 1L, prob = transition_matrix(model, a, s, sparse = FALSE))
        r <- reward_matrix(model, a, s, s_prime)
        a_prime <- greedy_action(Q, s_prime, epsilon)

        if (verbose) {
          if (i == 1L) {
            cat("\n*** Episode", e, "***\n")
          }
          cat(sprintf("Ep %i step %i - s:%s a:%i r:%.2f s':%s a':%i\n", 
                      e, i, s, a, r, s_prime, a_prime))
          #print(Q)
        }

        if (method == "sarsa") {
          # is called Sarsa because it uses the sequence s, a, r, s', a'
          Q[s, a] <-
            Q[s, a] + alpha * (r + gamma * Q[s_prime, a_prime] - Q[s, a])
          if (is.na(Q[s, a])) {
            Q[s, a] <- -Inf
          }
        } else if (method == "q_learning") {
          # a' is greedy instead of using the behavior policy
          a_max <- greedy_action(Q, s_prime, epsilon = 0)
          Q[s, a] <-
            Q[s, a] + alpha * (r + gamma * Q[s_prime, a_max] - Q[s, a])
          if (is.na(Q[s, a])) {
            Q[s, a] <- -Inf
          }
        } else if (method == "expected_sarsa") {
          p <-
            greedy_action(Q, s_prime, epsilon, prob = TRUE)
          exp_Q_prime <-
            sum(greedy_action(Q, s_prime, epsilon, prob = TRUE) * Q[s_prime, ],
              na.rm = TRUE
            )
          Q[s, a] <-
            Q[s, a] + alpha * (r + gamma * exp_Q_prime - Q[s, a])
          if (is.na(Q[s, a])) {
            Q[s, a] <- -Inf
          }
        }

        s <- s_prime
        a <- a_prime

        if (absorbing_states(model, state = s))
          break

        if (i >= horizon)
          break

        i <- i + 1L
      }
    }
    
    if (progress)
      pb$terminate()

    on.exit()
    
    model$solution <- list(
      method = method,
      alpha = alpha,
      epsilon = epsilon,
      n = n,
      Q = Q,
      converged = NA,
      policy = list(greedy_policy(Q))
    )

    model
  }
