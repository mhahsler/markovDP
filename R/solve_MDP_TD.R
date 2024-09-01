# Solve MDPs using Temporal Differencing

#' @rdname solve_MDP
#' @param alpha step size in `(0, 1]`.
#' @param epsilon used for \eqn{\epsilon}-greedy policies.
#' @param N number of episodes used for learning.
#' @export
solve_MDP_TD <-
  function(model,
           method = "q_learning",
           horizon = NULL,
           discount = NULL,
           alpha = 0.5,
           epsilon = 0.1,
           N = 100,
           U = NULL,
           progress = TRUE,
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
    gamma <- discount
    model$discount <- discount

    S <- model$states
    S_absorbing <- S[which(absorbing_states(model))]
    A <- model$actions
    P <- transition_matrix(model, sparse = TRUE)
    start <- .translate_belief(NULL, model = model, sparse = FALSE)

    method <-
      match.arg(method, c("sarsa", "q_learning", "expected_sarsa"))

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
      pb <- my_progress_bar(N)
    
    # loop episodes
    for (e in seq(N)) {
      if (progress)
        pb$tick()
      
      s <- sample(S, 1, prob = start)
      a <- greedy_action(Q, s, epsilon)

      # loop steps in episode
      i <- 1L
      while (TRUE) {
        s_prime <- sample(S, 1L, prob = P[[a]][s, ])
        r <- reward_matrix(model, a, s, s_prime)
        a_prime <- greedy_action(Q, s_prime, epsilon)

        if (verbose) {
          if (i == 1L) {
            cat("\n*** Episode", e, "***\n")
          }
          cat("Step", i, "- s a r s' a':", s, a, r, s_prime, a_prime, "\n")
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

        if (s %in% S_absorbing) {
          break
        }

        if (i >= horizon) {
          if (warn_horizon) {
            warning(
              "Max episode length of ",
              i,
              " reached without reaching an absorbing state! Terminating episode."
            )
          }
          break
        }

        i <- i + 1L
      }
    }

    model$solution <- list(
      method = method,
      alpha = alpha,
      epsilon = epsilon,
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
