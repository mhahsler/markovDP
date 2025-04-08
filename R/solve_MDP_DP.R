#' Solve MDPs using Dynamic Programming
#'
#' Solve MDPs via policy and value iteration.
#'
#' @family solver
#' 
#' @details
#' The following dynamic programming methods 
#' are implemented using the algorithms presented in Russell and Norvig (2010). 
#'
#' * **(Modified) Policy Iteration** (Howard 1960; Puterman and Shin 1978)
#' starts with a random policy and iteratively performs
#' a sequence of
#'   1. (Approximate) policy evaluation to estimate the value function for the
#'      current policy. Iterative policy evaluation can be approximated by 
#'      stopping early after `k_backups` iterations 
#'      (see [`policy_evaluation()`]. In this case the algorithm is called 
#'      _modified_ policy iteration.
#'   2. Policy improvement is performed by updating the policy to be greedy 
#'      (see [`greedy_policy()`])
#'      with respect to the new value function.
#' The algorithm stops when it converges to a stable policy (i.e., no changes
#' between two iterations). Note that the policy typically stabilizes before 
#' the value function converges. 
#'
#' * **Value Iteration** (Bellman 1957) starts with
#'   an arbitrary value function (by default all 0s) and iteratively
#'   updates the value function for each state using the Bellman update 
#'   equation (see [`bellman_update()`]).
#'   
#'   \deqn{v(s) \leftarrow \max_{a \in \mathcal{A}(s)} \sum_{s'} p(s' | s,a) [r(s,a, s') + \gamma v(s')]}
#'   
#'   The iteration
#'   is terminated when the solution converges or the maximum of `n` iterations
#'   has been reached.
#'   Approximate convergence is achieved
#'   for discounted problems (with \eqn{\gamma < 1})
#'   when the maximal value function change for any state \eqn{\delta} is
#'   \deqn{\delta \le \frac{error (1-\gamma)}{\gamma}.}
#'   It can be shown that this means
#'   that no state value is more than
#'   \eqn{error} from the value in the optimal value function. For undiscounted
#'   problems, we use \eqn{\delta \le error}.
#'
#'   A greedy policy
#'   is extracted from the final value function. Value iteration can be seen as
#'   policy iteration with policy evaluation truncated to one step.
#'
#' * **Prioritized Sweeping** (Moore and Atkeson, 1993; Andre et al., 1997; Li and Littman, 2008)
#'   approximate the optimal value
#'   function by iteratively adjusting one state at a time. While value and policy iteration
#'   sweep in every iteration through all states, prioritized sweeping 
#'   updates states in the order given by their priority.
#'   The priority reflects how much a state value may change
#'   given the most recently updated other states that can be directly reached via an action.
#'   This update order often lead to faster convergence compared
#'   to sweeping the whole state space in regular value iteration.
#'
#'   We implement the two priority update strategies described as __PS__ and
#'   __GenPS__ by Li and Littman (2008).
#'
#'   * __PS__ (Moore and Atkeson, 1993) updates the priority of a state \eqn{H(s)}
#'      using:
#'      \deqn{
#'        \forall{s \in \mathcal{S}}: H_{t+1}(s)  \leftarrow \begin{cases}
#'          \max(H_{t}(s), \Delta_t \max_{a \in \mathcal{A}}(p(s_t|s,a)) \text{ for } s \ne s_{t+1} \\
#'          \Delta_t \max_{a \in A}(p(s_t|s,a) \text{ for } s = s_{t+1}
#'          \end{cases}
#'      }
#'
#'      where \eqn{\Delta_t = |V_{t+1}(s_t) - V_t(s_t)| = |E(s_t; V_{t+1})|}, i.e.,
#'      the Bellman error for the updated state.
#'
#'   * __GenPS__ (Andre et al., 1997) updates all state priorities using their
#'      current Bellman error:
#'
#'      \deqn{\forall{s \in \mathcal{S}}: H_{t+1}(s) \leftarrow |E(s; V_{t+1})|}
#'
#'      where \eqn{E(s; V_{t+1}) = \max_{a \in A} \left[R(s,a) + \gamma \sum_{s \in S} p(s'|s,a) V(s')\right] - V(s)}
#'      is a state's Bellman error.
#'
#'   The update method can be chosen using the additional parameter `H_update`
#'   as the character string `"PS_random"`, `"PS_error"` or `"GenPS"`.
#'   The default is `H_update = "GenPS"`. For PS, random means that the
#'   priority vector is initialized with random values (larger than 0),
#'   and error means they are initialized with the Bellman error as in
#'   GenPS. However, this requires one complete sweep over all states.
#'
#'   This implementation stops updating when the largest priority values
#'   over all states is less than the specified `error`.
#'
#'   Since the algorithm does not sweep through the whole state space for each
#'   iteration, `n` is converted into an equivalent number of state updates
#'   \eqn{n = n\ |S|}.
#' 
#' @references
#' Andre, D., Friedman, N., and Parr, R. 1997. "Generalized prioritized sweeping." In Advances in Neural Information Processing Systems 10, pp. 1001-1007. [NeurIPS Proceedings](https://proceedings.neurips.cc/paper_files/paper/1997/file/7b5b23f4aadf9513306bcd59afb6e4c9-Paper.pdf)
#' 
#' Bellman, Richard. 1957. "A Markovian Decision Process." Indiana University Mathematics Journal 6: 679-84. [https://www.jstor.org/stable/24900506](https://www.jstor.org/stable/24900506).
#' 
#' Howard, R. A. 1960. "Dynamic Programming and Markov Processes." Cambridge, MA: MIT Press.
#' 
#' Li, Lihong, and Michael Littman. 2008. "Prioritized Sweeping Converges to the Optimal Value Function." DCS-TR-631. Rutgers University. \doi{10.7282/T3TX3JSX}
#' 
#' Moore, Andrew, and C. G. Atkeson. 1993. "Prioritized Sweeping: Reinforcement Learning with Less Data and Less Real Time." Machine Learning 13 (1): 103–30. \doi{10.1007/BF00993104}.
#' 
#' Puterman, Martin L., and Moon Chirl Shin. 1978. "Modified Policy Iteration Algorithms for Discounted Markov Decision Problems." Management Science 24: 1127-37. \doi{10.1287/mnsc.24.11.1127}.
#' 
#' Russell, Stuart J., and Peter Norvig. 2020. Artificial Intelligence: A Modern Approach (4th Edition). Pearson. [http://aima.cs.berkeley.edu/](http://aima.cs.berkeley.edu/).
#' 
#' @examples
#' data(Maze)
#' 
#' maze_solved <- solve_MDP(Maze, method = "DP:VI", verbose = TRUE)
#' policy(maze_solved)
#' 
#' # use prioritized sweeping (which is known to be fast for mazes)
#' maze_solved <- solve_MDP(Maze, method = "DP:GenPS", verbose = TRUE)
#' policy(maze_solved)
#'
#' # finite horizon
#' maze_solved <- solve_MDP(Maze, method = "DP:VI", horizon = 3)
#' policy(maze_solved)
#' gw_plot(maze_solved, epoch = 1)
#' gw_plot(maze_solved, epoch = 2)
#' gw_plot(maze_solved, epoch = 3)
#' 
#' @inheritParams solve_MDP
#' @param method string; one of the following solution methods:
#'    * `'VI'` - value iteration
#'    * `'PI'` - policy iteration
#'    * `'GenPS'`, `'PS_error'`, `'PS_random'` - prioritized sweeping 
#' @param V a vector with initial state values. If
#'   `NULL`, then the default of a vector of all 0s ([V_zero()]) is used.
#' @param n maximum number of iterations allowed to converge. If the
#'   maximum is reached then the non-converged solution is returned with a
#'   warning.
#' @param error value iteration: maximum Bellman error allowed for the 
#'    convergence criterion.
#' @param k_backups policy iteration: maximum number of Bellman backups used in 
#'    the iterative policy evaluation step. Policy evaluation typically converges earlier
#'    with a maximum Bellman error less than `error`.  
#' 
#' @inherit solve_MDP return 
#' 
#' @export
solve_MDP_DP <- function(model,
                         method = "VI",
                         horizon = NULL,
                         discount = NULL,
                         n = 1000L,
                         error = 0.001,
                         k_backups = 10L,
                         V = NULL,
                         ...,
                         matrix = TRUE,
                         continue = FALSE,
                         verbose = FALSE,
                         progress = TRUE) {
  .nodots(...)
  
  methods <- c("VI",
               "PI",
               "GenPS",
               "PS_error",
               "PS_random")
  method <- match.arg(method, methods)

  if (!inherits(model, "MDP"))
    stop("This model requires transition probabilities in the MDP description.")
  
  model <- .prep_model(model, horizon, discount, matrix, verbose, progress)
   
  if (continue) {
    if (is.null(model$solution$policy[[1]]$V))
      stop("model solution does not contain a V vector to continue from!")
    V <- model$solution$policy[[1]]$V
  }
  
  ret <- switch(
    method,
    VI = {
      if (is.infinite(model$horizon)) {
        solve_MDP_DP_VI(
          model,
          error,
          n,
          V = V,
          progress = progress,
          verbose = verbose,
          ...
        )
      } else {
        solve_MDP_DP_VI_finite_horizon(
          model,
          horizon = model$horizon,
          V = V,
          progress = progress,
          verbose = verbose,
          ...
        )
      }
    },
    PI = {
      if (is.infinite(model$horizon)) {
        solve_MDP_DP_PI(
          model,
          n,
          k_backups,
          V = V,
          progress = progress,
          verbose = verbose,
          ...
        )
      } else {
        stop("Method not implemented yet for finite horizon problems.")
      }
    },
    PS_error = {
      if (is.infinite(model$horizon)) {
        solve_MDP_PD_PS(
          model,
          error = error,
          n = n,
          V = V,
          H_update = "PS_error",
          progress = progress,
          verbose = verbose,
          ...
        )
      } else {
        stop("Method not implemented yet for finite horizon problems.")
      }
    },
    PS_random = {
      if (is.infinite(model$horizon)) {
        solve_MDP_PD_PS(
          model,
          error = error,
          n = n,
          V = V,
          H_update = "PS_random",
          progress = progress,
          verbose = verbose,
          ...
        )
      } else {
        stop("Method not implemented yet for finite horizon problems.")
      }
    },
    GenPS = {
      if (is.infinite(model$horizon)) {
        solve_MDP_PD_PS(
          model,
          error = error,
          n = n,
          V = V,
          H_update = "GenPS",
          progress = progress,
          verbose = verbose,
          ...
        )
      } else {
        stop("Method not implemented yet for finite horizon problems.")
      }
    }
  )
  
  ret$solution$method <- method
  
  ret
}

# TODO: we could check for convergence
solve_MDP_DP_VI_finite_horizon <-
  function(model,
           horizon,
           V = NULL,
           progress = TRUE,
           verbose = FALSE,
           ...) {
    .nodots(...)
    
    if (!inherits(model, "MDP")) {
      stop("'model' needs to be of class 'MDP'.")
    }
    
    if (progress) {
      pb <- my_progress_bar(horizon, name = "solve_MDP")
      pb$tick(0)
    }
    
    horizon <- as.integer(horizon)
    policy <- vector(mode = "list", length = horizon)
    
    if (is.null(V)) 
      V <- V_zero(model)
    
    if (verbose) {
      cat("Running value iteration (finite horizon)")
      cat("\nhorizon: ", horizon)
      cat("\nInitial V (first 20 max):\n")
      print(head(V, n = 20))
      
      cat("\n")
    }
    
    for (t in seq(horizon, 1)) {
      if (progress)
        pb$tick()
      
      if (verbose > 1) {
        cat("Iteration for t = ", t, "\n")
      }
      
      V_prime <- bellman_update(model, V)
      pi <- V_prime$pi
      V <- V_prime$V
      
      policy[[t]] <- .normalize_policy(pi, model, V = V)
    }
    
    if (progress)
      pb$terminate()
    
    model$solution <- list(policy = policy, converged = NA)
    model
  }

solve_MDP_DP_VI <- function(model,
           error,
           n = 1000L,
           V = NULL,
           progress = TRUE,
           verbose = FALSE,
           ...) {
    .nodots(...)
    
    if (verbose)
      progress <- FALSE
    
    if (progress) {
      pb <- my_progress_spinner(name = "solve_MDP", format_extra = " | max delta U: :delta | press esc/CTRL-C to terminate early")
      pb$tick(0, token = list(delta = "-"))
    }
    
    discount <- model$discount
    if (discount < 1) {
      convergence_factor <- (1 - discount) / discount
    } else {
      convergence_factor <- 1
    }
    convergence_limit <- convergence_factor * error
    
    if (is.null(V)) 
      V <- V_zero(model)
    
    # return unconverged result when interrupted
    on.exit({
      if (progress) {
        pb$tick(0, token =
                  list(delta = paste0(
                    signif(delta, 3), "/", signif(convergence_limit, 3)
                  )))
        pb$terminate()
      }
      
      if (!converged) {
        if (i < n)
          warning("Manual interupt: MDP solver did not converge at iteration ", i)
        else
          warning(
            "MDP solver did not converge after ",
            i,
            " iterations.",
            " Consider decreasing the 'discount' factor or increasing 'error' or 'n'."
          )
      } 
      
      if (verbose) {
        cat("Terminated at iteration:", i, "\n")
      }
      
      model$solution <- list(
        method = "value iteration",
        policy = list(.normalize_policy(pi, model, V = V)),
        converged = converged,
        delta = delta,
        iterations = i
      )
      return(model)
    })
   
    if (verbose) {
      cat("Running value iteration (infinite horizon)")
      cat("\nerror for convergence:  ", error) 
      cat("\nn (max iterations):     ", n) 
      cat("\nInitial V (first 20 max):\n")
      print(head(V, n = 20))
      
      cat("\n")
    } 
    
    converged <- FALSE
    delta <- Inf
    for (i in seq_len(n)) {
      if (progress)
        pb$tick(token =
                  list(delta = paste0(
                    signif(delta, 3), "/", signif(convergence_limit, 3)
                  )))
      
      if (verbose > 1) {
        cat("Iteration:", i)
      }
      
      V_prime <- bellman_update(model, V)
      pi <- V_prime$pi
      V_prime <- V_prime$V
      
      delta <- max(abs(V_prime - V))
      V <- V_prime
      
      if (verbose > 1) {
        cat("\t -> delta:", delta, "\n")
      }
      
      if (delta <= convergence_limit) {
        converged <- TRUE
        break
      }
    }
    
    # return via on.exit()
  }

## Policy iteration
solve_MDP_DP_PI <-
  function(model,
           n = 1000,
           k_backups = 10,
           V = NULL,
           progress = TRUE,
           verbose = FALSE,
           ...) {
    .nodots(...)
    
    if (progress) {
      pb <- my_progress_spinner(
        name = "solve_MDP",
        format_extra = paste0(
          " | changed actions: :changed_actions/",
          length(S(model)),
          " | press esc/CTRL-C to terminate early"
        )
      )
      pb$tick(0, token = list(changed_actions = "-"))
    }
    
    S <- S(model)
    
    if (is.null(V))
      V <- V_zero(model)
    
    pi <- greedy_policy(Q_values(model, V))
    V <- pi$V
    pi <- pi$action
    
    if (progress)
      pb$tick(0, token = list(changed_actions = "-"))
    
    # return unconverged result when interrupted
    on.exit({
      if (progress) {
        pb$tick(0, token = list(changed_actions = changed_actions))
        pb$terminate()
      }
      
      if (!converged) {
        if (i < n)
          warning("Manual interupt: MDP solver did not converge at iteration ", i)
        else
          warning(
            "MDP solver did not converge after ",
            i,
            " iterations.",
            " Consider decreasing the 'discount' factor or increasing 'n'."
          )
      } 
      
      if (verbose) {
        cat("\nTerminated at iteration:", i, "\n")
      }
      
      model$solution <- list(
        method = "policy iteration",
        policy = list(.normalize_policy(pi, model, V = V)),
        converged = converged,
        iterations = i
      )
      
      return(model)
    })
    
    if (verbose) {
      cat("Running policy iteration (infinite horizon)")
      cat("\nn (max updates):       ", n) 
      cat("\nk_backups:             ", k_backups) 
      cat("\nInitial V (first 20 max):\n")
      print(head(V, n = 20))
      
      cat("\n")
    } 
    
    changed_actions <- length(S(model))
    converged <- FALSE
    i <- 0L
    while (i < n) {
      i <- i + 1L
      if (progress)
        pb$tick(token = list(changed_actions = changed_actions))
      
      # evaluate to get V from pi
      V <- policy_evaluation(
        model,
        pi,
        V,
        k_backups = k_backups,
        progress = FALSE
      )
      
      V_prime <- bellman_update(model, V)
      pi_prime <- V_prime$pi
      V <- V_prime$V
      Q <- V_prime$Q
      
      # changed_actions <- sum(pi != pi_prime)
      # account for random tie breaking using Q
      changed_actions <- sum(Q[cbind(seq_len(nrow(Q)), pi)] != apply(Q, MARGIN = 1, max))
      
      if (verbose > 1)
        cat(i, "| # of updated actions ", changed_actions, "\n")
      
      if (changed_actions == 0) {
        converged <- TRUE
        break
      }
      
      pi <- pi_prime
    }
    
    # return is done in on.exit()
  }

## Prioritized sweeping
## Details are in: L. Li and M.L. Littman: Prioritized sweeping
##    converges to the optimal value function. Technical report DCSTR-631,
##    Department of Computer Science, Rutgers University, May 2008
## https://www.academia.edu/15223291/Prioritized_Sweeping_Converges_to_the_Optimal_Value_Function
##
## We initialize H(s) using the states reward so we start with the biggest
## reward states and propagate the reward backwards.
solve_MDP_PD_PS <- function(model,
           error,
           n = 1000L,
           V = NULL,
           H_update = "GenPS",
           use_n_times_states = TRUE,
           progress = TRUE,
           verbose = FALSE,
           ...) {
    .nodots(...)
    
    S <- S(model)
    A <- A(model)
    discount <- model$discount
    
    if (use_n_times_states)
      n <- as.integer(n * length(S))
    
    if (is.na(n))
      stop("n * |S| is >", .Machine$integer.max, "! Reduce n.")
     
    if (verbose > 1) 
      progress <- FALSE
    
    i <- 0L
    
    if (progress) {
      pb <- my_progress_spinner(name = "solve_MDP",
                                ticks = "state updates",
                                format_extra = " | max priority: :error | press esc/CTRL-C to terminate early")
      pb$tick(0, token = list(error = "-"))
    }
    
    H_update <- match.arg(H_update, c("GenPS", "PS_random", "PS_error"))
    
    H <- NULL
    Q <- NULL
    
    # this may be a continue if we have H, Q and pi from the previous run
    if (!is.null(V)) {
      H <- model$solution$H
      Q <- model$solution$Q
      pi <- as.integer(policy(model)$action)
    } else {
      V <- rep(0, times = length(S))
    }
    
    # needs initialization
    if (is.null(H)) {
      if (H_update == "PS_random") {
        pi <- as.integer(random_policy(model)$action)
        H <- runif(length(S)) + 1e-6 + error
        Q <- NULL
      } else {   ### GenPS or PS_error
        i <- i + length(S)
        VQpi <- bellman_update(model, V)
        pi <- as.integer(VQpi$pi)
        
        # Note: PS may not find the optimal value if some priorities are 0 (See paper)
        H <- abs(VQpi$V - V) + 1e-6
        V <- VQpi$V
        
        if (H_update == "GenPS")
          Q <- VQpi$Q
        else
          Q <- NULL
      }
    } else {
      if (verbose)
        cat("Continuing...\n")
    }
    
    # absorbing states get 0 priority
    #H[absorbing_states(model, sparse = "states")] <- 0
    
    # return result also when interrupted
    on.exit({
      if (progress) {
        pb$tick(0, tokens = list(error = paste0(signif(err, 3), "/", signif(error, 3))))
        pb$terminate()
      }
         
      if (!converged) {
        if (i < n)
          warning("Manual interupt: MDP solver did not converge at state update ", i)
        else
          warning(
            "MDP solver did not converge after ",
            i,
            " iterations (delta = ",
            signif(delta, 3),
            ").",
            " Consider decreasing the 'discount' factor or increasing 'error' or 'n'."
          )
      } 
       
      if (verbose) {
        cat("State updates performed:", i, "(equivalent to", signif(i/length(S), 3), 
            "complete sweeps)", "\n")
        cat("Converged:", converged, "\n")
      }
      
      model$solution <- list(
        method = "prioritized sweeping",
        policy = list(.normalize_policy(pi, model, V = V)),
        Q = Q,
        H = H,
        converged = converged,
        state_updates = i
      )
      return(model)
    })
    
    if (verbose) {
      cat("Running prioritized sweeping")
      cat("\nH_update:              ", H_update) 
      cat("\nerror for convergence: ", error) 
      cat("\nn (max updates):       ", n) 
      cat("\nInitial V (first 20 max):\n")
      print(head(V, n = 20))
      
      cat("\n")
    } 
    
    err <- sum(H)
    delta <- Inf
    converged <- FALSE
    while (i < n) {
      i <- i + 1L
      if (progress) {
        pb$tick(tokens = list(error = paste0(signif(err, 3), "/", signif(error, 3))))
      }
        
      # update state with highest priority (error) next
      s_t <- which.max.random(H)
      
      if (verbose > 1)
        cat(i, "  updating state", s_t, sQuote(S[s_t]), "| action", pi[s_t], "-> ")
      
      ### do Bellman update for a single state
      V_prime <- .bellman_state_update(model, s_t, V)
      
      pi[s_t] <- V_prime$pi
      V_prime <- V_prime$V
      delta <- V_prime - V[s_t]
      
      if (verbose > 1)
        cat(pi[s_t], "| V ", signif(V[s_t], 3), "->", V_prime, "\n")
      
      V[s_t] <- V_prime
      
      
      # Bellman error in state s: E(s; V) =
      #      max_a [r(s,a) + gamma sum_{s in S} p(s'|s,a) V(s')] - V(s)
      
      if (H_update != "GenPS") {
        # This is Moore and Atkeson (1993) called PS in Li and Littman (2008)
        #
        # forall s: H_t+1(s) <- max(H_t(s), |delta| max_a(p(s_t|s,a)) for s != s_t+1
        #                       delta max_a(p(s_t|s,a))) for s = s_t+1
        #
        # where delta = V_t+1(s_t) - V_t(s_t) = E(s_t; V_t+1)
        #
        # Note: changes only happen if p(s_t|s,a) > 0 for any a
        
        H[s_t] <- 0 # it will be updated with delta max... below
        
        # find states that can get us to s, i.e., max_a((P(s_t|s, a)) > 0
        max_a_T <- apply(
          transition_matrix(
            model,
            end.state = s_t,
            simplify = TRUE,
            sparse = FALSE
          ),
          MARGIN = 1,
          max
        )
        
        s_update <- which(max_a_T > 0)
        
        # Note: this will also take care of s_t
        H[s_update] <- pmax(H[s_update], abs(delta) * max_a_T[s_update])
      }
      
      else {  #if (H_update == "GenPS") {
        # This is Andre, Friedman, & Parr (1998) called GenPS in Li and Littman (2008)
        #
        # for all s: H_t+1 <- |E(s; V_t+1)|
        #
        # with E(s; V_t+1) = max_a [r(s,a) + gamma sum_{s in S} r(s'|s,a) V(s')] - V(s)
        #
        # for s != s_t: 
        #    E(s; V_t+1) = max_a [E^a(s; V_t) + gamma p(s_t|s,a) E(s_t; V_t)]     
        # for s = s_t:
        #    E(s_t; V_t+1) = max_a [E^a(s_t; V_t) + gamma p(s_t|s_t,a) E(s_t; V_t)] - E(s_t; V_t)
        #
        # Notes:
        #  * |E(s_t; V_t)| = H(s_t)_t
        #  * changes only happen for states s where p(s_t|s,a) > 0 for any a
        # We need: E^a(s; V_t) for all these states and E(s_t; V_t)
        # E^a(s; .) is Q

        p_sa <- transition_matrix(
          model,
          end.state = s_t,
          simplify = TRUE,
          sparse = FALSE
        )
        
        s_update <- union(which(rowSums(p_sa) > 0), s_t)
        E_st <- max(Q[s_t, ])
        
        Q[s_update, ] <- Q[s_update, ] + discount * p_sa[s_update, ] * E_st
        Q[s_t, ] <- Q[s_t, ] - E_st
        
        H[s_update] <- abs(apply(Q[s_update, , drop = FALSE], MARGIN = 1, max))
      }
      
      err <- max(H)
      
      if (verbose > 1) {
        cat(
          "    updating H for states:",
          paste(s_update, sQuote(S[s_update]), "->", signif(H[s_update], 3), collapse = ", "),
          "\n"
        )
        cat("    max priority (error):", err, ">=", error, "\n")
      }
      
      if (err < error) {
        converged <- TRUE
        break
      }
    }
  
    # return via on.exit()  
  }


