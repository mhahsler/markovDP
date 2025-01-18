#' Solve MDPs using Temporal Differencing
#'
#' Solve MDPs using 1-step and n-step tabular temporal difference control 
#' methods like q-learning and Sarsa.
#'
#' @family solver
#' 
#' @details
#' Implemented are several tabular temporal difference control methods
#' described in Sutton and Barto (2018).
#' Note that the MDP transition and reward models are used
#' for these reinforcement learning methods only to sample from
#' the environment.
#' 
#' The implementation uses an \eqn{\epsilon}-greedy behavior policy,
#' where the parameter `epsilon` controls the degree of exploration. 
#' The algorithms use a step size parameter \eqn{\alpha} (learning rate). 
#' The learning rate `alpha` can be specified as a 
#' function with the signature `function(t, n)`, where `t` is the number of episodes 
#' processed and `n` is the number of updates for the entry in the Q-table. 
#' 
#' The general 1-step update is
#' \deqn{
#' Q(S_t,A_t) \leftarrow Q(S_t,A_t) + \alpha [G_t - Q(S_t,A_t)],
#' }
#' where \eqn{G_t} is the target estimate for the given q-value. The different methods below
#' estimate the target value differently.
#'
#' If the model has absorbing states to terminate episodes, then no maximal episode length
#' (`horizon`) needs to
#' be specified. To make sure that the algorithm does finish in a reasonable amount of time,
#' episodes are stopped after 1000 actions (with a warning). For models without absorbing states,
#' the episode length has to be specified via `horizon`.
#'
#' * **Q-Learning** (Watkins and Dayan 1992) is an off-policy temporal difference method that uses
#'    an \eqn{\epsilon}-greedy behavior policy and learns a greedy target
#'    policy. The target value is estimated as the one-step bootstrapping using the 
#'    target greedy policy: 
#'    \deqn{G_t = R_{t+1} + \gamma \max_a Q(S_{t+1}, a)}
#'
#' * **Sarsa** (Rummery and Niranjan 1994) is an on-policy method that follows and learns
#'    the same policy. Here a an \eqn{\epsilon}-greedy policy is used. 
#'    The final \eqn{\epsilon}-greedy policy is converted into a greedy policy.
#'    \eqn{\epsilon} can be lowered over time (see parameter `continue`) 
#'    to learn a greedy policy. The target is estimated 
#'    as the one-step bootstrapping following the behavior policy:
#'    \deqn{G_t = R_{t+1} + \gamma Q(S_{t+1}, A_{t+1})}
#'    
#'
#' * **Expected Sarsa** (Sutton and Barto 2018). We implement an on-policy 
#'   Sarsa with an \eqn{\epsilon}-greedy policy which uses the
#'   the expected value under the current policy for the update.
#'   It moves deterministically in the same direction as Sarsa would
#'   move in expectation. 
#'   \deqn{G_t = R_{t+1} + \gamma \sum_a \pi(a|S_{t+1})Q(S_{t+1}, a)}
#'   
#'   
#'   Because it uses the expectation, we can
#'   set the step size \eqn{\alpha} to large values and 1 is common.
#'   The off-policy use of expected Sarsa simplifies to 
#'   the Q-learning algorithm.  
#'   
#' * **On and off-policy n-step Sarsa** (Sutton and Barto 2018).
#'   Estimate the return using the last \eqn{n} time steps:
#'   \deqn{
#'   G_{t:t+n} = R_{t+1} + \gamma R_{t+2} + ... + \gamma^{n-1} R_{t+n} + \gamma^n Q(S_{t+n}, A_{t+n})
#'   }
#'   
#'   \eqn{n = 1} is regular 1-step Sarsa, \eqn{n = \inf} is equivalent to 
#'   Monte Carlo Control.
#'   
#'   This estimate is used as the target for Sarsa. For the off-policy case,
#'   the update uses the importance sampling ratio. Note that updates are delayed 
#'   \eqn{n} steps in this backward looking algorithm.
#' 
#' @references 
#'  
#' Rummery, G., and Mahesan Niranjan. 1994. "On-Line Q-Learning Using Connectionist Systems." Techreport CUED/F-INFENG/TR 166. Cambridge University Engineering Department.
#' 
#' Sutton, R. 1988. "Learning to Predict by the Method of Temporal Differences." Machine Learning 3: 9-44. [https://link.springer.com/article/10.1007/BF00115009](https://link.springer.com/article/10.1007/BF00115009).
#'
#' Sutton, Richard S., and Andrew G. Barto. 2018. Reinforcement Learning: An Introduction. Second. The MIT Press. [http://incompleteideas.net/book/the-book-2nd.html](http://incompleteideas.net/book/the-book-2nd.html).
#' 
#' Watkins, Christopher J. C. H., and Peter Dayan. 1992. "Q-Learning." Machine Learning 8 (3): 279-92. \doi{10.1007/BF00992698}.
#' 
#' @inheritParams solve_MDP
#' @param alpha step size as a function of the time step `t` and the number of times
#'   the respective Q-value was updated `n` or a scalar. For expected Sarsa, alpha is
#'   often set to 1.
#' @param epsilon used for the \eqn{\epsilon}-greedy behavior policies.
#' @param n number of episodes used for learning.
#' @param Q an initial state-action value matrix. By default an all 0 matrix is 
#'        used.
#' 
#' @inherit solve_MDP return 
#' 
#' @examples
#' data(Maze)
#' 
#' # Example 1: Learn a Policy using Q-Learning
#' maze_learned <- solve_MDP(Maze, method = "TD:q_learning",
#'     epsilon = 0.2, n = 500, horizon = 100, verbose = TRUE)
#' maze_learned
#'
#' policy(maze_learned)
#' plot_value_function(maze_learned)
#' gw_plot(maze_learned)
#'
#' # Keep on learning, but with a reduced epsilon
#' maze_learned <- solve_MDP(maze_learned, method = "TD:q_learning",
#'     epsilon = 0.01, n = 500, horizon = 100, continue = TRUE, verbose = TRUE)
#'
#' policy(maze_learned)
#' plot_value_function(maze_learned)
#' gw_plot(maze_learned)
#'
#' # Example 2: n-step Sarsa
#' maze_learned <- solve_MDP(Maze, method = "TDN:sarsa_on_policy",
#'     n_step = 3, n = 10, horizon = 100, verbose = TRUE)
#' maze_learned
#' 
#' gw_plot(maze_learned)
#' @export
solve_MDP_TD <-
  function(model,
           method = "q_learning",
           horizon = NULL,
           discount = NULL,
           alpha = function(t, n)
             min(10 / n, 1),
           epsilon = 0.2,
           n = 1000,
           Q = 0,
           ...,
           matrix = TRUE,
           continue = FALSE,
           progress = TRUE,
           verbose = FALSE) {
    .nodots(...)
    
    method <-
      match.arg(method, c("sarsa", "q_learning", "expected_sarsa"))
    
    model <- .prep_model(model, horizon, discount, matrix, verbose, progress)
     
    if (!is.finite(model$horizon)) {
      stop("Finite horizon needed for exploration!")
    }
    
    if (!is.function(alpha))
      alpha_val <- alpha
    
    if (progress) {
      pb <- my_progress_bar(n + 1L, name = "solve_MDP")
      pb$tick(0)
    }
    
    # Initialize Q and Q_N
    if (continue) {
      Q <- model$solution$Q
      Q_N <- model$solution$Q_N
      if (is.null(Q) || is.null(Q_N))
        stop("model solution does not contain a Q matrix to continue from!")
    } else {
      Q <- init_Q(model, Q)
      Q_N <-
        matrix(0L,
               nrow = nrow(Q),
               ncol = ncol(Q),
               dimnames = dimnames(Q))
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
    
    S <- S(model)
    A <- A(model)
    start <- start_vector(model, sparse = FALSE)
    horizon <- model$horizon
    discount <- model$discount
     
    # loop episodes
    e <- 0L
    while (e < n) {
      e <- e + 1L
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
        
        # act
        a_res <- act.int(model, s, a)
        s_prime <- a_res$state_prime 
        r <- a_res$r 
        
        # act without a function call
        #s_prime <- sample(S, 1L, prob = transition_matrix(model, a, s, sparse = FALSE))
        #r <- reward_matrix(model, a, s, s_prime)
       
        
        # for Sarsa we need (s, a, r, s', a')
        a_prime <- greedy_action(Q, s_prime, epsilon)
        
        if (verbose > 1) {
          if (i == 1L) {
            cat("\n*** Episode", e, "***\n")
          }
          cat(
            sprintf(
              "Ep %i step %i - s=%s a=%i r=%.3f s'=%s a'=%i - q(s,a): %.3f -> ",
              e,
              i,
              s,
              a,
              r,
              s_prime,
              a_prime,
              Q[s, a]
            )
          )
        }
        
        # update Q and Q_N and calculate alpha
        Q_N[s, a] <- Q_N[s, a] + 1L
        if (is.function(alpha))
          #alpha_val <- alpha(t, Q_N[s, a])
          alpha_val <- alpha(t, sum(Q_N[s, ]))
        
        target <- switch(
          method,
          # on-policy: used s' and a' from the behavior policy
          sarsa = Q[s_prime, a_prime],
          
          # off-policy: uses an estimate of the the target greedy policy (max(Q))
          q_learning = max(Q[s_prime, ]),
          
          # on-policy: Uses expectation under the behavior policy
          # (off-policy would use the expectation under the greedy behavior policy -> q-learning)
          expected_sarsa = sum(greedy_action(Q, s_prime, epsilon, prob = TRUE) * Q[s_prime, ], na.rm = TRUE)
        )
        
        Q[s, a] <- Q[s, a] + alpha_val * (r + discount * target - Q[s, a])
        
        if (is.na(Q[s, a])) {
          Q[s, a] <- -Inf
        }
        
        if (verbose > 1) {
          cat(sprintf("%.3f (N: %i alpha: %.3f)\n", Q[s, a], Q_N[s, a], alpha_val))
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


#' @rdname solve_MDP_TD
#' @param n_step integer; steps for bootstrapping
#' @export
solve_MDP_TDN <-
  function(model,
           method = "sarsa_on_policy",
           horizon = NULL,
           discount = NULL,
           n_step,
           alpha = function(t, n)
             min(10 / n, 1),
           epsilon = 0.2,
           n = 1000,
           Q = 0,
           matrix = TRUE,
           continue = FALSE,
           progress = TRUE,
           verbose = FALSE,
           ...) {
    .nodots(...)
    
    if (missing(n_step))
      stop("argument \"n_step\" is missing!")
    
    method <-
      match.arg(method,
                c("sarsa_on_policy", "sarsa_off_policy"))
   
    model <- .prep_model(model, horizon, discount, matrix, verbose, progress)
    
    if (!is.finite(model$horizon)) {
      stop("Finite horizon needed for exploration!")
    }
    
    if (!is.function(alpha))
      alpha_val <- alpha
    
    if (progress) {
      pb <- my_progress_bar(n + 1L, name = "solve_MDP")
      pb$tick(0)
    }
    
    
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
               dimnames = dimnames(Q))
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
      cat("\nn_step            ", n_step, "\n")
      
      cat("\nInitial Q (first 20 max):\n")
      print(head(Q, n = 20))
      
      cat("\nInitial Q_N (first 20 max):\n")
      print(head(Q_N, n = 20))
      cat("\n")
    }
    
    S <- S(model)
    A <- A(model)
    start <- start_vector(model, sparse = FALSE)
    horizon <- model$horizon
    discount <- model$discount
    
    # store S_t, A_t and R_t; initialize as 0; accessed with t %% n_step + 1
    S_t <- integer(n_step + 1L)
    A_t <- integer(n_step + 1L)
    R_t <- numeric(n_step + 1L)
    gamma_t <- numeric(n_step + 1L)
    gamma_n_steps <- discount^n_step
    
    # index in circular buffer
    t2idx <- function(t) t %% (n_step + 1) + 1
   
    # loop episodes
    e <- 0L
    while (e < n) {
      e <- e + 1L
      if (progress)
        pb$tick()
      
      s <- sample.int(length(S), 1L, prob = start)
      a <- greedy_action(Q, s, epsilon)
      
      S_t[1L] <- s
      A_t[1L] <- a
      
      # loop steps in episode t = 0, 1, ... tt
      t <- 0L
      tt <- Inf
      discount_t <- 1
      while (TRUE) {
        t_plus_1_idx <- t2idx(t + 1)
        if ((t %% 100 == 0) && progress && !pb$finished)
          pb$tick(0)
        
        if (t < tt) {
          # take action a
          a_res <- act.int(model, s, a)
          s_prime <- a_res$state_prime 
          r <- a_res$r 
          
          # no funciton call
          #s_prime <- sample(S, 1L, prob = transition_matrix(model, a, s, sparse = FALSE))
          #r <- reward_matrix(model, a, s, s_prime)
          
          S_t[t_plus_1_idx] <- s_prime
          R_t[t_plus_1_idx] <- r
          gamma_t[t_plus_1_idx] <- discount_t <- discount_t * discount # discount^t
          
          if (absorbing_states(model, state = s_prime)) {
            tt <- t + 1   # end of episode T
          } else {
            # choose an action using the behavior policy (which is epsilon greedy here)
            a_prime <- greedy_action(Q, s_prime, epsilon)
            A_t[t_plus_1_idx] <- a_prime
          }
        } else {
          # lean up memory
          s_prime <- S_t[t_plus_1_idx] <- NA_integer_
          r <- R_t[t_plus_1_idx] <- 0
        }
        
        tau <- t - n_step + 1 # time step that is updated (this is n_steps back)
        
        if (verbose > 1) {
          if (t == 0L) {
            cat("\n*** Episode", e, "***\n")
          }
          cat(
            sprintf(
              "Step %i: tau=%i s=%s a=%i r=%.3f s'=%s a'=%i",
              t,
              tau,
              s,
              a,
              r,
              s_prime,
              a_prime
            )
          )
        }
        
        if (tau >= 0) {
          G <- sum(gamma_t * R_t)
          
          # add final reward estimate
          if (tau + n_step < tt) {
            tau_plus_n_idx <- t2idx(tau + n_step)
            G <- G + gamma_n_steps * Q[S_t[tau_plus_n_idx], A_t[tau_plus_n_idx]]
          }
          
          tau_idx <-  t2idx(tau)
          s_tau <- S_t[tau_idx]
          a_tau <- A_t[tau_idx]
          
          # update Q and Q_N and calculate alpha
          Q_N[s_tau, a_tau] <- Q_N[s_tau, a_tau] + 1L
          if (is.function(alpha))
            #alpha_val <- alpha(t, Q_N[s_tau, a_tau])
            alpha_val <- alpha(t, sum(Q_N[s_tau, ]))
          
          if (verbose > 1) {
            cat(sprintf("; q(%s,%i):%.3f ->", s_tau, a_tau, Q[s_tau, a_tau]))
          }
          
          # on-policy: use 1
          rho <- 1
          
          if (method == "n_step_SARSA_off_policy") {
            # off-policy: use importance sampling ratio
            # rho_t:t+n-1 = prod_k=t^min(t+n-1,T) pi(A_k|S_k)/b(A_k|S_k)
            # only updates if all actions would have been chosen by pi!
            pi_hits <- (apply(Q, MARGIN = 1, which.max.random)[S_t] == A_t)
            pi_hits <- pi_hits[!is.na(pi_hits)]
            if (all(pi_hits, na.rm = TRUE)) {
              # greedy pi / epsilon-greedy b
              rho <- 1 / (1 - epsilon + epsilon / length(S))^length(pi_hits)
            } else {
              rho <- 0
            }
          }
          
          Q[s_tau, a_tau] <- Q[s_tau, a_tau] + alpha_val * rho * (G - Q[s_tau, a_tau])
          
          if (verbose > 1) {
            cat(sprintf(" %.3f (N: %i alpha: %.3f, rho: %.3f)", Q[s_tau, a_tau], Q_N[s_tau, a_tau], alpha_val, rho))
          }
          
        }
          
        if (verbose > 1) {
          cat("\n")
          }
          
        if (tau >= tt - 2)
          break
        
        s <- s_prime
        a <- a_prime
        
        t <- t + 1L
      }
    }
    
    # return via on.exit()
  }