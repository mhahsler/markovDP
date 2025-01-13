# Solve MDPs using Dynamic Programming (Value and policy iteration)

#' @rdname solve_MDP
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
#' @param continue logical; Continue with an unconverged solution specified in `model`.
#' @param matrix logical; if `TRUE` then matrices for the transition model and
#'    the reward function are taken from the model first. This can be slow if functions
#'    need to be converted or do not fit into memory if the models are large. If these
#'    components are already matrices, then this is very fast. For `FALSE`, the
#'    transition probabilities and the reward is extracted when needed. This is slower,
#'    but removes the time and memory requirements needed to calculate the matrices.
#' @export
solve_MDP_DP <- function(model,
                         method = "value_iteration",
                         horizon = NULL,
                         discount = NULL,
                         n = 1000L,
                         error = 0.001,
                         k_backups = 10L,
                         V = NULL,
                         matrix = TRUE,
                         continue = FALSE,
                         progress = TRUE,
                         verbose = FALSE,
                         ...) {
  if (!inherits(model, "MDP")) {
    stop("'model' needs to be of class 'MDP'.")
  }
  
  if (verbose > 1)
    progress <- FALSE

  methods <- c("value_iteration",
               "policy_iteration",
               "prioritized_sweeping",
               "GenPS",
               "PS_error",
               "PS_random")
  method <- match.arg(method, methods)
  
  ### default is infinite horizon
  model$horizon <- horizon %||% model$horizon %||% Inf
  
  model$discount <- discount %||% model$discount
  if (is.null(model$discount)) {
    message("No discount rate specified. Using .9 for the infinite horizon problem!")
    model$discount <- .9
  }
  
  if (continue) {
    if (is.null(model$solution$policy[[1]]$V))
      stop("model solution does not contain a V vector to continue from!")
    V <- model$solution$policy[[1]]$V
  }
  
  if (matrix) {
    if (verbose)
      cat("Precomputing matrices for R and T ...")
    model <- normalize_MDP(
      model,
      sparse = NULL,
      precompute_absorbing = FALSE,
      progress = progress
    )
    
    if (verbose)
      cat(" done.\n")
  }
  
  ret <- switch(
    method,
    value_iteration = {
      if (is.infinite(model$horizon)) {
        MDP_value_iteration_inf_horizon(
          model,
          error,
          n,
          V = V,
          progress = progress,
          verbose = verbose,
          ...
        )
      } else {
        MDP_value_iteration_finite_horizon(
          model,
          horizon = model$horizon,
          V = V,
          progress = progress,
          verbose = verbose,
          ...
        )
      }
    },
    policy_iteration = {
      if (is.infinite(model$horizon)) {
        MDP_policy_iteration_inf_horizon(
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
    prioritized_sweeping = {
      if (is.infinite(model$horizon)) {
        MDP_PS_inf_horizon(
          model,
          error = error,
          n = n,
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
        MDP_PS_inf_horizon(
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
        MDP_PS_inf_horizon(
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
        MDP_PS_inf_horizon(
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
MDP_value_iteration_finite_horizon <-
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

MDP_value_iteration_inf_horizon <-
  function(model,
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

## Prioritized sweeping
## Details are in: L. Li and M.L. Littman: Prioritized sweeping
##    converges to the optimal value function. Technical report DCSTR-631,
##    Department of Computer Science, Rutgers University, May 2008
## https://www.academia.edu/15223291/Prioritized_Sweeping_Converges_to_the_Optimal_Value_Function
##
## We initialize H(s) using the states reward so we start with the biggest
## reward states and propagate the reward backwards.
MDP_PS_inf_horizon <-
  function(model,
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


## Policy iteration

MDP_policy_iteration_inf_horizon <-
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
