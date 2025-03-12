#' Solve MDPs using Tabular Temporal Differencing
#'
#' Solve MDPs using tabular temporal difference control
#' methods like q-learning and 1-step, n-step Sarsa, and Sarsa(\eqn{\lambda}).
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
#' The `epsilon` and the learning rate `alpha` can be specified as 
#' fixed numbers between 0 and 1 or as a schedule function where the 
#' value is gradually reduced (see Schedule section below).
#'
#' All temporal differencing methods us the following update: 
#' \deqn{
#' Q(S_t,A_t) \leftarrow Q(S_t,A_t) + \alpha [G_t - Q(S_t,A_t)],
#' }
#' where \eqn{G_t} is the target estimate for the given q-value. The different methods below
#' estimate the target value differently (see 1 and n-step methods below).
#' 
#' If the model has absorbing states to terminate episodes, then no maximal episode length
#' (`horizon`) needs to
#' be specified. To make sure that the algorithm does finish in a reasonable amount of time,
#' episodes are stopped after 1000 actions (with a warning). For models without absorbing states,
#' the episode length has to be specified via `horizon`.
#'
#' 
#' ## 1-step Methods
#' 
#' 1-step methods estimate the return by just looking one state ahead.
#' 
#' * **Q-Learning** (Watkins and Dayan 1992) is an off-policy (learns a greedy target policy) 
#'    temporal difference method. We use here an \eqn{\epsilon}-greedy behavior policy. 
#'    The update target value is estimated by one-step bootstrapping using the
#'    reward and the value of the state following the current greedy target policy:
#'    \deqn{G_t = R_{t+1} + \gamma \max_a Q(S_{t+1}, a)}
#'
#' * **Sarsa** (Rummery and Niranjan 1994) is an on-policy method (behavior and 
#'    target policy are the same). We use an \eqn{\epsilon}-greedy policy
#'    and the final \eqn{\epsilon}-greedy policy is converted into a greedy policy.
#'    \eqn{\epsilon} can be lowered over time (see [schedule] and parameter `continue`)
#'    to learn a approximately greedy policy. The target is estimated
#'    as the one-step bootstrap estimate following the current behavior policy:
#'    \deqn{G_t = R_{t+1} + \gamma Q(S_{t+1}, A_{t+1})}
#'
#' * **Expected Sarsa** (Sutton and Barto 2018) learns the behavior policy 
#'   (on-policy learning).
#'   We use Sarsa with an \eqn{\epsilon}-greedy policy which uses the
#'   the expected value under the current policy for the update:
#'   \deqn{G_t = R_{t+1} + \gamma \sum_a \pi(a|S_{t+1})Q(S_{t+1}, a)}
#'   
#'   Expected Sarsa moves in the same direction as Sarsa would
#'   move in expectation.
#'   Because it uses the expectation, we can
#'   set the step size \eqn{\alpha} to large values and 1 is common.
#'   
#'   The Q-learning algorithm can be seen as a simplification of the 
#'   off-policy version of expected Sarsa.
#'
#' ## n-step Methods
#'
#' n-step methods use a longer look ahead. The 
#' return is estimated by looking 
#' \eqn{n} time steps ahead (`n_step` in the code) and using 
#' the rewards and then the value of the reached state:
#'   
#' \deqn{
#'   G_{t:t+n} = R_{t+1} + \gamma R_{t+2} + ... + \gamma^{n-1} R_{t+n} + \gamma^n Q(S_{t+n}, A_{t+n})
#' }
#'   
#' While n-step methods conceptually look ahead, the implementation has to wait 
#' for the values to become available. This means that updates are \eqn{n} 
#' steps delayed, i.e., the update for step \eqn{t} is performed at \eqn{t+n}.
#'   
#' * **n-step Sarsa** (Sutton and Barto 2018).
#'   The estimated return is used as the update target for Sarsa. 
#'   
#'   `n_step = 1` is regular 1-step Sarsa. Using `n_step = Inf` is equivalent to
#'   Monte Carlo Control, however, [`solve_MDP_MC()`] is more memory efficient.  
#'
#'   Sarsa learns on-policy (i.e., the behavioral policy which is often 
#'   an exploring \eqn{\epsilon}-greedy policy). For the off-policy case,
#'   when the optimal policy is learned, the update are corrected using the 
#'   importance sampling ratio. 
#'
#' ## Eligibility Traces 
#' 
#' 
#' Eligibility traces also look ahead but the update is performed in a backward looking manner
#' by using a memory parameter called the trace that remember what states or Q-values lead to the 
#' current event. We focus here on tabular Sarsa(\eqn{\lambda})
#' 
#' From the forward-view perspective, the \eqn{\lambda}-reward can be seen as 
#' an average of infinite n-step backups:
#' \deqn{G_t^\lambda = (1-\lambda) \sum_{n=1}^\infty \lambda^{n-1} R_t^{(n)},}
#' where \eqn{\lambda} is the eligibility trace decay factor. The implementation 
#' uses the backward view with the TD-error
#' \deqn{\delta_t = R_{t+1} + \gamma Q(S_{t+1}, A_{t+1}) - Q(S_{t}, A_{t})}
#' and the update
#' \deqn{Q(s,a) = Q(s,a) + \alpha \delta_t z_t(s,a) \qquad \forall s, a.}  
#' 
#' The trace is updated at each time step as
#' \deqn{z_t(s,a) = \gamma z_{t-1}(s,a) \qquad \forall s, a}
#' and \eqn{z_t(S_t,A_t) = z_t(S_t,A_t) + 1}. This means, the trace is set to one
#' for the current state and action and then decays with the parameter \eqn{\lambda}.
#' 
#'
#' ## Schedules
#' 
#' * epsilon schedule: `t` is increased by each processed episode.
#' * alpha schedule: `t` is set to the number of times the a Q-value for state 
#'      `s` was updated. 
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
#' @param method string; one of the following solution methods:
#'  `"sarsa"`, `"q_learning"`, or `"expected_sarsa"`
#' @param alpha step size (learning rate). A scalar value between 0 and 1 or a 
#'    [schedule].
#' @param epsilon used for the \eqn{\epsilon}-greedy behavior policies. A scalar value between 0 and 1 or a 
#'    [schedule].
#' @param n_step number of steps to look ahead for n-step Sarsa.
#' @param lambda eligibility trace decay factor for Sarsa(lambda).
#' @param on_policy logical; should we learn on-policy (the behavior policy) 
#'        vs. off-policy (using importance sampling)? Only used for method `"sarsa"`.
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
#'     n = 200, horizon = 100)
#' maze_learned
#'
#' policy(maze_learned)
#' gw_plot(maze_learned)
#'
#' # Example 2: Learn a Policy using 1-step Sarsa
#' maze_learned <- solve_MDP(Maze, method = "TD:sarsa",
#'     n = 200, horizon = 100)
#' maze_learned
#'
#' policy(maze_learned)
#' gw_plot(maze_learned)
#'
#' # Example 3: Perform one episode for 3-step Sarsa
#' 
#' # run one episode in verbose mode.
#' maze_learned <- solve_MDP(Maze, method = "TD:sarsa",
#'     n_step = 5, n = 1, horizon = 100, verbose = 2)
#'     
#' # verbose output:
#' #  * tau ... time step updated (laggs n_step)
#' #  * -> ... update of the Q-value
#' #  * N ... number of updates for this state
#' #  * G ... reward estimate for the n_steps
#' #  * alpha ... learning rate (schedule may depend on N)
#' #  * rho ... importance sampling ratio (1 for on-policy learning)
#'
#' # run more episode
#' maze_learned <- solve_MDP(Maze, method = "TD:sarsa",
#'     n_step = 5, n = 100, horizon = 100)
#' maze_learned
#' policy(maze_learned)
#' gw_plot(maze_learned)
#' 
#' # Example: Tabular Sarsa(lambda)
#' maze_learned <- solve_MDP(Maze, method = "TD:sarsa",
#'     lambda = .1, n = 1, horizon = 100)
#' maze_learned
#' policy(maze_learned)
#' gw_plot(maze_learned)
#' @export
solve_MDP_TD <-
  function(model,
           method = "sarsa",
           horizon = NULL,
           discount = NULL,
           alpha = schedule_exp(0.2, .01),
           epsilon = schedule_exp(1, 0.1),
           n_step = 1,
           lambda = 0,
           on_policy = TRUE,
           n,
           Q = 0,
           ...,
           matrix = TRUE,
           continue = FALSE,
           progress = TRUE,
           verbose = FALSE) {
    
    if (n_step == 1 && lambda == 0) {
      solve_MDP_TD_1_step(model,
                          method,
                          horizon,
                          discount,
                          alpha,
                          epsilon,
                          on_policy,
                          n,
                          Q,
                          ...,
                          matrix = matrix,
                          continue = continue,
                          progress = progress,
                          verbose = verbose)
    } else if (lambda == 0) {
      solve_MDP_TD_n_step(model,
                          method,
                          horizon,
                          discount,
                          alpha,
                          epsilon,
                          n_step,
                          on_policy,
                          n,
                          Q,
                          ...,
                          matrix = matrix,
                          continue = continue,
                          progress = progress,
                          verbose = verbose)
  } else if (lambda != 0) {
    if (n_step != 1)
      stop("you can only specify n_step or lambda, but not both!")
      solve_MDP_TD_lambda(model,
                          method,
                          horizon,
                          discount,
                          alpha,
                          epsilon,
                          lambda,
                          on_policy,
                          n,
                          Q,
                          ...,
                          matrix = matrix,
                          continue = continue,
                          progress = progress,
                          verbose = verbose)
    
    
    }
}
 
   
solve_MDP_TD_1_step <-
  function(model,
           method = "sarsa",
           horizon = NULL,
           discount = NULL,
           alpha = schedule_exp(0.2, .01),
           epsilon = schedule_exp(1, 0.1),
           on_policy,
           n,
           Q = 0,
           ...,
           matrix = TRUE,
           continue = FALSE,
           progress = TRUE,
           verbose = FALSE) {
    .nodots(...)
    
    method <-
      match.arg(method, c("q_learning", "sarsa", "expected_sarsa"))
    
    # this works for MDP and MDPTF with a defined state space
    if (!inherits(model, "MDPE") || is.null(S(model)))
      stop("The model needs to be an MDP description with a specified state space.")
    
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
    
    model <- .prep_model(model, horizon, discount, matrix, verbose, progress)
    
    if (!is.finite(model$horizon)) {
      stop("Finite horizon needed for exploration!")
    }
    
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
        n_step = 1L,
        on_policy = on_policy,
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
      cat("\nalpha:            ", show_schedule(alpha))
      cat("\nepsilon:          ", show_schedule(epsilon))
      cat("\nn                 ", n, "\n")
    }  
     
    if (verbose > 1) {   
      cat("\nInitial Q (first 20 max):\n")
      print(head(Q, n = 20))
      
      cat("\nInitial Q_N (first 20 max):\n")
      print(head(Q_N, n = 20))
      cat("\n")
    }
    
    S <- S(model)
    A <- A(model)
    horizon <- model$horizon
    gamma <- model$discount
    
    # loop episodes
    e <- 0L
    while (e < n) {
      e <- e + 1L
      if (progress)
        pb$tick()
      
      ### epsilon func
      if (!is.null(epsilon_seq)) 
        epsilon <- epsilon_seq[e]
      ###
      
      # get episode start state and first action
      s <- start(model, as = "id")
      a <- greedy_action(model, s, Q, epsilon)
      
      # loop steps in episode
      i <- 0L
      while (TRUE) {
        i <- i + 1L
        if ((i %% 100 == 0) && progress && !pb$finished)
          pb$tick(0)
        
        # act
        a_res <- act(model, s, a, fast = TRUE)
        s_prime <- a_res$state_prime
        r <- a_res$r
        
        # MDPTF: features -> id
        if (is.matrix(s_prime)) s_prime <- normalize_state_id(s_prime, model)
      
        # for Sarsa we need (s, a, r, s', a')
        a_prime <- greedy_action(model, s_prime, Q, epsilon)
        
        if (verbose > 1) {
          if (i == 1L) {
            cat("\n*** Episode", e, "***\n")
          }
          cat(
            sprintf(
              "Ep %i step %i\ts=%-10s a=%-10s r=%.3f\ts'=%-10s a'=%-10s\tq(s,a): %.3f -> ",
              e,
              i,
              normalize_state_label(s, model),
              normalize_action_label(a, model),
              r,
              normalize_state_label(s_prime, model),
              normalize_action_label(a_prime, model),
              Q[s, a]
            )
          )
        }
        
        # update Q and Q_N and calculate alpha
        Q_N[s, a] <- Q_N[s, a] + 1L
        
        ### alpha/epsilon func
        # alpha uses the count and not the episode number!
        if (!is.null(alpha_func)) 
          alpha <- alpha_func(sum(Q_N[s, ]))
        ###
        
        G <- switch(
          method,
          
          # off-policy: uses an estimate of the the target greedy policy (max(Q))
          q_learning = max(Q[s_prime, ]),
          
          # on-policy: used s' and a' from the behavior policy
          sarsa = Q[s_prime, a_prime],
          
          # on-policy: Uses expectation under the behavior policy
          # (off-policy would use the expectation under the greedy behavior policy -> q-learning)
          expected_sarsa = sum(greedy_action(model, s_prime, Q, epsilon, prob = TRUE) * Q[s_prime, ], na.rm = TRUE)
        )
        
        # on-policy: always use importance sampling factor of 1
        rho <- 1
        
        ## off_policy for sarsa
        if (!on_policy && method != "q_learning")  {
          # off-policy: use importance sampling ratio
          # rho = prob greedy pi / prob epsilon-greedy b
          if (as.integer(a_prime) %in% which(Q[s_prime, ] == max(Q[s_prime, ]))) {
            rho <-  1 / (1 - epsilon + epsilon / length(S))
          } else {
            rho <- 0
          }
        }
        
        Q[s, a] <- Q[s, a] + alpha * rho * (r + gamma * G - Q[s, a])
        
        if (is.na(Q[s, a])) {
          Q[s, a] <- -Inf
        }
        
        if (verbose > 1) {
          cat(sprintf("%.3f (N: %i alpha: %.3f)\n", Q[s, a], Q_N[s, a], alpha))
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


solve_MDP_TD_n_step <-
  function(model,
           method = "sarsa_on_policy",
           horizon = NULL,
           discount = NULL,
           alpha = schedule_exp(0.2, .01),
           epsilon = schedule_exp(1, 0.1),
           n_step,
           on_policy = TRUE,
           n,
           Q = 0,
           ...,
           matrix = TRUE,
           continue = FALSE,
           progress = TRUE,
           verbose = FALSE
           ) {
    .nodots(...)
    
    method <-
      match.arg(method, c("sarsa"))
    
    if (missing(n_step))
      stop("argument \"n_step\" is missing!")
    
    model <- .prep_model(model, horizon, discount, matrix, verbose, progress)
    
    if (!is.finite(model$horizon)) {
      stop("Finite horizon needs to be specified to avoid potential infinite loops!")
    }
    
    
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
    
    # n_step = Inf: this is MC control. We size all the arrays to fit horizon
    if (!is.finite(n_step)) {
      n_step = horizon + 1
    }
    
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
        on_policy = on_policy,
        n_step = n_step,
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
      cat("\nalpha:            ", show_schedule(alpha))
      cat("\nepsilon:          ", show_schedule(epsilon))
      cat("\nn                 ", n, "\n")
      cat("\nn_step            ", n_step, "\n")
    }
    
    if (verbose > 1) {
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
    gamma <- model$discount
    
    # store S_t, A_t and R_t; initialize as 0; accessed with t %% n_step + 1
    S_t <- integer(n_step + 1L)
    A_t <- integer(n_step + 1L)
    R_t <- numeric(n_step + 1L)
    gamma_t <- numeric(n_step + 1L)
    gamma_n_steps <- gamma^n_step
    
    # index in circular buffer
    t2idx <- function(t)
      t %% (n_step + 1) + 1
    
    # loop episodes
    e <- 0L
    while (e < n) {
      e <- e + 1L
      if (progress)
        pb$tick()
      
      ### epsilon func
      if (!is.null(epsilon_seq)) 
        epsilon <- epsilon_seq[e]
      ###
      
      s <- sample.int(length(S), 1L, prob = start)
      a <- greedy_action(model, s, Q, epsilon)
      
      S_t[1L] <- s
      A_t[1L] <- a
      
      # loop steps in episode t = 0, 1, ... tt
      t <- 0L
      tt <- Inf
      gamma_tt <- 1
      while (TRUE) {
        t_plus_1_idx <- t2idx(t + 1)
        if ((t %% 100 == 0) && progress && !pb$finished)
          pb$tick(0)
        
        if (t < tt) {
          # take action a
          a_res <- act(model, s, a, fast = TRUE)
          s_prime <- a_res$state_prime
          r <- a_res$r
          
          # no function call version (may be faster)
          #s_prime <- sample(S, 1L, prob = transition_matrix(model, a, s, sparse = FALSE))
          #r <- reward_matrix(model, a, s, s_prime)
          
          S_t[t_plus_1_idx] <- s_prime
          R_t[t_plus_1_idx] <- r
          gamma_t[t_plus_1_idx] <- gamma_tt <- gamma_tt * gamma # gamma^t
          
          if (absorbing_states(model, state = s_prime)) {
            tt <- t + 1   # end of episode T
          } else {
            # choose an action using the behavior policy (which is epsilon greedy here)
            a_prime <- greedy_action(model, s_prime, Q, epsilon)
            A_t[t_plus_1_idx] <- a_prime
          }
        } else {
          # clean up memory
          #s_prime <- S_t[t_plus_1_idx] <- NA_integer_
          r <- R_t[t_plus_1_idx] <- 0
        }
        
        tau <- t - n_step + 1L # time step that is updated (this is n_steps back)
        
        if (verbose > 1) {
          if (t == 0L) {
            cat("\n*** Episode", e, "***\n")
          }
          cat(
            sprintf(
              "Step %i (tau=%i)\ts=%-10s a=%-10s r=%.3f\ts'=%-10s a'=%-10s",
              t,
              tau,
              normalize_state_label(s, model),
              normalize_action_label(a, model),
              r,
              normalize_state_label(s_prime, model),
              normalize_action_label(a_prime, model)
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
          
          ### alpha/epsilon func
          # alpha uses the count and not the episode number!
          if (!is.null(alpha_func)) 
            alpha <- alpha_func(sum(Q_N[s, ]))
          ###
          
          if (verbose > 1) {
            cat(sprintf("-> q(%s,%s):%.3f ->", normalize_state_label(s_tau, model), 
                        normalize_action_label(a_tau, model), Q[s_tau, a_tau]))
          }
          
          # on-policy: use 1
          rho <- 1
          
          if (!on_policy) {
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
          
          Q[s_tau, a_tau] <- Q[s_tau, a_tau] + alpha * rho * (G - Q[s_tau, a_tau])
          
          if(is.na(Q)[s_tau, a_tau]) {
            #browser()
            stop("A NA Q value has been generated! Please report this with a code example as a bug.")
          }
          
          
          if (verbose > 1) {
            cat(sprintf(
              " %.3f (N=%i G=%.3f alpha=%.3f rho=%.2f)",
              Q[s_tau, a_tau],
              Q_N[s_tau, a_tau],
              G,
              alpha,
              rho
            ))
          }
          
        }
        
        if (verbose > 1) {
          cat("\n")
        }
        
        if (tau >= tt)
          break
        
        s <- s_prime
        a <- a_prime
        
        t <- t + 1L
      }
    }
    
    # return via on.exit()
  }


solve_MDP_TD_lambda <-
  function(model,
           method = "sarsa",
           horizon = NULL,
           discount = NULL,
           alpha = schedule_exp(0.2, .01),
           epsilon = schedule_exp(1, 0.1),
           lambda, 
           on_policy,
           n,
           Q = 0,
           ...,
           matrix = TRUE,
           continue = FALSE,
           progress = TRUE,
           verbose = FALSE) {
    .nodots(...)
    
    method <-
      match.arg(method, c("sarsa"))
    
    # this works for MDP and MDPTF with a defined state space
    if (!inherits(model, "MDPE") || is.null(S(model)))
      stop("The model needs to be an MDP description with a specified state space.")
    
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
    
    model <- .prep_model(model, horizon, discount, matrix, verbose, progress)
    
    if (!is.finite(model$horizon)) {
      stop("Finite horizon needed for exploration!")
    }
    
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
    
    # eligibility vector (as an s x a matrix)
    z <- matrix(0L,
                nrow = nrow(Q),
                ncol = ncol(Q),
                dimnames = dimnames(Q))
    
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
        lambda = lambda,
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
      cat("\nalpha:            ", show_schedule(alpha))
      cat("\nepsilon:          ", show_schedule(epsilon))
      cat("\nlambda:          ", lambda)
      cat("\nn                 ", n, "\n")
    }  
    
    if (verbose > 1) {   
      cat("\nInitial Q (first 20 max):\n")
      print(head(Q, n = 20))
      
      cat("\nInitial Q_N (first 20 max):\n")
      print(head(Q_N, n = 20))
      cat("\n")
    }
    
    S <- S(model)
    A <- A(model)
    horizon <- model$horizon
    gamma <- model$discount
    
    # loop episodes
    e <- 0L
    while (e < n) {
      e <- e + 1L
      if (progress)
        pb$tick()
      
      ### epsilon func
      if (!is.null(epsilon_seq)) 
        epsilon <- epsilon_seq[e]
      ###
      
      # get episode start state and first action
      s <- start(model, as = "id")
      a <- greedy_action(model, s, Q, epsilon)
      
      # loop steps in episode
      i <- 0L
      while (TRUE) {
        i <- i + 1L
        if ((i %% 100 == 0) && progress && !pb$finished)
          pb$tick(0)
        
        # act
        a_res <- act(model, s, a, fast = TRUE)
        s_prime <- a_res$state_prime
        r <- a_res$r
        
        # MDPTF: features -> id
        if (is.matrix(s_prime)) s_prime <- normalize_state_id(s_prime, model)
        
        # for Sarsa we need (s, a, r, s', a')
        a_prime <- greedy_action(model, s_prime, Q, epsilon)
        
        if (verbose > 1) {
          if (i == 1L) {
            cat("\n*** Episode", e, "***\n")
          }
          cat(
            sprintf(
              "Ep %i step %i\ts=%-10s a=%-10s r=%.3f\ts'=%-10s a'=%-10s\tq(s,a): %.3f ->",
              e,
              i,
              normalize_state_label(s, model),
              normalize_action_label(a, model),
              r,
              normalize_state_label(s_prime, model),
              normalize_action_label(a_prime, model),
              Q[s, a]
            )
          )
        }
        
        # update Q and Q_N and calculate alpha
        Q_N[s, a] <- Q_N[s, a] + 1L
        
        ### alpha/epsilon func
        # alpha uses the count and not the episode number!
        if (!is.null(alpha_func)) 
          alpha <- alpha_func(sum(Q_N[s, ]))
        ###
        
        delta <- r + gamma * Q[s_prime, a_prime] - Q[s,a]
        
        z[s,a] <- z[s,a] + 1 
        
        # FIXME: maybe we need individual alphas!
        Q <- Q + alpha * delta * z
         
        ## off_policy for sarsa
        if (!on_policy)  {
          # importance sampling ratio
          stop("Off-policy Sarsa not implemented yet.")
        }
        
        if (verbose > 1) {
          cat ( sprintf(
            " %.3f\n",
            Q[s, a]
          ))
          cat("z:\n")
          print(z)
        }
        
        # decay aligibility trace vector
        z <- gamma * lambda * z 
        
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




