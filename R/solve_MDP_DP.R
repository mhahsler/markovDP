# Solve MDPs using Dynamic Programming (Value and policy iteration)

#' @rdname solve_MDP
#' @param U a vector with initial utilities used for each state. If
#'   `NULL`, then the default of a vector of all 0s is used.
#' @param iter_max maximum number of iterations allowed to converge. If the
#'   maximum is reached then the non-converged solution is returned with a
#'   warning.
#' @param error value iteration: maximum error allowed in the utility of any state
#'   (i.e., the maximum policy loss) used as the termination criterion.
#' @param k_backups policy iteration: number of look ahead steps used for approximate policy evaluation
#'   used by the policy iteration method.
#' @param continue logical; Continue with an unconverged solution specified in `model`.
#' @param matrix logical; if `TRUE` then matrices for the transition model and
#'    the reward function are taken from the model first. This can be slow if functions
#'    need to be converted or not fit into memory if the models are large. If these
#'    components are already matrices, then this is very fast. For `FALSE`, the
#'    transition probabilities and the reward is extracted when needed. This is slower,
#'    but removes the time to calculate the matrices and it saves memory.
#' @export
solve_MDP_DP <- function(model,
                         method = "value_iteration",
                         horizon = NULL,
                         discount = NULL,
                         iter_max = 1000,
                         error = 0.01,
                         k_backups = 10,
                         U = NULL,
                         matrix = TRUE,
                         continue = FALSE,
                         progress = TRUE,
                         verbose = FALSE,
                         ...) {
  if (!inherits(model, "MDP")) {
    stop("'model' needs to be of class 'MDP'.")
  }
  
  if (verbose)
    progress <- FALSE
  
  # Note: This is used by gridworld_animate()
  dots <- list(...)
  if (!is.null(dots$n))
    iter_max <- dots$n
  
  methods <- c("value_iteration",
               "policy_iteration",
               "prioritized_sweeping")
  method <- match.arg(method, methods)
  
  ### default is infinite horizon
  if (!is.null(horizon)) {
    model$horizon <- horizon
  }
  if (is.null(model$horizon)) {
    model$horizon <- Inf
  }
  
  if (!is.null(discount)) {
    model$discount <- discount
  }
  if (is.null(model$discount)) {
    message("No discount rate specified. Using .9 for the infinite horizon problem!")
    model$discount <- .9
  }
  
  if (continue) {
    if (is.null(model$solution$policy[[1]]$U))
      stop("model solution does not contain a U vector to continue from!")
    U <- model$solution$policy[[1]]$U
  }
  
  if (matrix) {
    if (verbose)
      cat("Precomputing dense matrices for R and T ...")
    model <- normalize_MDP(
      model,
      sparse = FALSE,
      precompute_absorbing = FALSE,
      precompute_unreachable = FALSE,
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
          iter_max,
          U = U,
          progress = progress,
          verbose = verbose,
          ...
        )
      } else {
        MDP_value_iteration_finite_horizon(
          model,
          horizon = model$horizon,
          U = U,
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
          iter_max,
          k_backups,
          U = U,
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
          iter_max = iter_max,
          U = U,
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
           U = NULL,
           progress = TRUE,
           verbose = FALSE,
           ...) {
    if (!inherits(model, "MDP")) {
      stop("'model' needs to be of class 'MDP'.")
    }
    
    if (progress) {
      pb <- my_progress_bar(horizon, name = "solve_MDP")
      pb$tick(0)
    }
    
    S <- model$states
    horizon <- as.integer(horizon)
    
    if (is.null(U)) {
      U <- rep(0, times = length(S))
    }
    names(U) <- S
    
    policy <- vector(mode = "list", length = horizon)
    
    for (t in seq(horizon, 1)) {
      if (progress)
        pb$tick()
      
      if (verbose) {
        cat("Iteration for t = ", t, "\n")
      }
      
      U_prime <- bellman_update(model, U)
      pi <- U_prime$pi
      U <- U_prime$U
      
      policy[[t]] <- data.frame(
        state = S,
        U = U,
        action = pi,
        row.names = NULL
      )
    }
    
    if (progress)
      pb$terminate()
    
    model$solution <- list(policy = policy, converged = NA)
    model
  }

MDP_value_iteration_inf_horizon <-
  function(model,
           error,
           N_max = 1000,
           U = NULL,
           matrix = TRUE,
           progress = TRUE,
           verbose = FALSE,
           ...) {
    if (verbose)
      progress <- FALSE
    
    if (progress) {
      pb <- my_progress_spinner(name = "solve_MDP", format_extra = " | max delta U: :delta | press esc/CTRL-C to terminate early")
      pb$tick(0, token = list(delta = "-"))
    }
    
    S <- model$states
    
    GAMMA <- model$discount
    if (GAMMA < 1) {
      convergence_factor <- (1 - GAMMA) / GAMMA
    } else {
      convergence_factor <- 1
    }
    convergence_limit <- convergence_factor * error
    
    
    if (is.null(U)) {
      U <- rep(0, times = length(S))
    }
    names(U) <- S
    
    # return unconverged result when interrupted
    on.exit({
      warning("MDP solver interrupted early.")
      
      if (verbose) {
        cat("\nTerminated during iteration:", i, "\n")
      }
      
      model$solution <- list(
        method = "value iteration",
        policy = list(data.frame(
          state = S,
          U = U,
          action = pi,
          row.names = NULL
        )),
        converged = converged,
        delta = delta,
        iterations = i
      )
      return(model)
    })
    
    converged <- FALSE
    delta <- Inf
    for (i in seq_len(N_max)) {
      if (progress)
        pb$tick(token =
                  list(delta = paste0(
                    signif(delta, 3), "/", signif(convergence_limit, 3)
                  )))
      
      if (verbose && !progress) {
        cat("Iteration:", i)
      }
      
      U_t_plus_1 <- bellman_update(model, U)
      pi <- U_t_plus_1$pi
      U_t_plus_1 <- U_t_plus_1$U
      
      delta <- max(abs(U_t_plus_1 - U))
      U <- U_t_plus_1
      
      if (verbose && !progress) {
        cat(" -> delta:", delta, "\n")
      }
      
      if (delta <= convergence_limit) {
        converged <- TRUE
        break
      }
    }
    
    if (progress)
      pb$tick(0, token =
                list(delta = paste0(
                  signif(delta, 3), "/", signif(convergence_limit, 3)
                )))
    
    if (verbose) {
      cat("\nIterations needed:", i, "\n")
      # print(U)
    }
    
    if (!converged) {
      warning(
        "MDP solver did not converge after ",
        i,
        " iterations (delta = ",
        delta,
        ").",
        " Consider decreasing the 'discount' factor or increasing 'error' or 'N_max'."
      )
    }
    
    # clear on.exit
    on.exit()
    
    model$solution <- list(
      method = "value iteration",
      policy = list(data.frame(
        state = S,
        U = U,
        action = pi,
        row.names = NULL
      )),
      converged = converged,
      delta = delta,
      iterations = i
    )
    model
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
           iter_max = 10000,
           U = NULL,
           H_init = "reward",
           H_update = "GenPS",
           matrix = TRUE,
           progress = TRUE,
           verbose = FALSE,
           ...) {
    if (verbose)
      progress <- FALSE
    
    if (progress) {
      pb <- my_progress_spinner(name = "solve_MDP",
                                ticks = "state updates",
                                format_extra = " | max priority: :error | press esc/CTRL-C to terminate early")
      pb$tick(0, token = list(error = "-"))
    }
    
    S <- model$states
    A <- model$actions
    GAMMA <- model$discount
    
    if (!is.null(U)) {
      H <- model$solution$H
      pi <- as.integer(model$solution$policy[[1]]$action)
      if (is.null(pi))
        stop("Policy missing in the model. Can't continue.")
      
    } else {
      U <- rep(0, times = length(S))
      pi <- as.integer(random_policy(model)$action)
      H <- NULL
    }
    names(U) <- S
    
    # initialize with the maximum reward. This is a single value iteration pass.
    H_init <- match.arg(H_init, c("reward", "random"))
    
    if (is.null(H)) {
      if (H_init == "reward") {
        if (verbose)
          cat("initializing priority H using the Belmann error.\n\n")
        H <- abs(bellman_update(model, U = U)$U) + error + 1e-6
      } else {
        ### random
        if (verbose)
          cat("initializing priority H randomly greater than error.\n\n")
        H <- runif(length(model$states)) + error + 1e-6
      }
      
      # cat(round(H,2), "\n")
    }
    # may not find the optimal value if some priorities are 0 (See paper)
    
    # absorbing states get 0 priority
    H[absorbing_states(model, sparse = "states")] <- 0
    
    
    # return unconverged result when interrupted
    on.exit({
      warning("MDP solver did not converge (manual interrupt).")
      
      if (verbose) {
        cat("\nTerminated during iteration:", i, "\n")
      }
      
      pi <- factor(pi, levels = seq_along(A), labels = A)
      
      model$solution <- list(
        method = "prioritized sweeping",
        policy = list(data.frame(
          state = S,
          U = U,
          action = pi,
          row.names = NULL
        )),
        converged = converged,
        delta = delta,
        H = H,
        iterations = i
      )
      return(model)
    })
    
    err <- sum(H)
    delta <- Inf
    converged <- FALSE
    for (i in seq_len(iter_max)) {
      if (progress)
        pb$tick(tokens = list(error = paste0(signif(err, 3), "/", signif(error, 3))))
      
      # update state with highest priority (error) next
      s_t <- which.max.random(H)
      
      if (verbose)
        cat(i, "  updating state", s_t, sQuote(S[s_t]), "| action", pi[s_t], "-> ")
      
      ### do Bellman update for a single state
      Us_prime <- .bellman_state_update(model, s_t, U)
      
      pi[s_t] <- Us_prime$pi
      Us_prime <- Us_prime$U
      delta <- abs(Us_prime - U[s_t])
      
      if (verbose)
        cat(pi[s_t], "| U ", signif(U[s_t], 3), "->", Us_prime, "\n")
      
      U[s_t] <- Us_prime
      
      
      # Bellman error in state s: E(s; V) =
      #      max_a [R(s,a) + gamma sum_{s in S} T(s'|s,a) V(s')] - V(s)
      
      if (H_update == "PS") {
        # This is Moore and Atkeson (1993) called PS in Li and Littman (2008)
        #
        # forall s: H_t+1(s) <- max(H_t(s), delta max_a(T(s_t|s,a))) for s != s_t+1
        #                       delta max_a(T(s_t|s,a))) for s = s_t+1
        #
        # where delta = |V_t+1(s_t) - V_t(s_t)| = |E(s; V_t+1)|
        
        H[s_t] <- 0 # it will be updated with delta max... below
        
        # find states that can get us to s, i.e., max_a((T(s_t|s, a)) > 0
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
        
        ancestor_ids <- which(max_a_T > 0)
        
        # note this also increases the priority for the state we just updated from
        H[ancestor_ids] <- pmax(H[ancestor_ids], delta * max_a_T[ancestor_ids])
      }
      
      else if (H_update == "GenPS") {
        # This is Andre, Friedman, & Parr (1998) called GenPS in Li and Littman (2008)
        #
        # for all s: H_t+1 <- |E(s; V_t+1)| =
        #             max_a [R(s,a) + gamma sum_{s in S} T(s'|s,a) V(s')] - V(s)
        #
        # Note: |E(s; V_t+1)| only changes for s_t and it's ancestors
        T_s_t <- transition_matrix(
          model,
          end.state = s_t,
          simplify = TRUE,
          sparse = FALSE
        )
        
        ancestor_ids <- union(which(rowSums(T_s_t) > 0), s_t)
  
        
        # note this also increases the priority for the state we just updated from
        H[ancestor_ids] <- apply((reward_matrix(model, start.state = ancestor_ids, end.state = s_t, simplify = TRUE) + 
                 GAMMA * T_s_t[ancestor_ids,] * U[s_t]) - U[ancestor_ids], MARGIN = 1, max)
        
      }
      
      err <- max(H)
      
      if (verbose) {
        cat(
          "    updating H for states:",
          paste(ancestor_ids, sQuote(S[ancestor_ids]), "->", H[ancestor_ids], collapse = ", "),
          "\n"
        )
        cat("    max priority (error):", err, ">=", error, "\n")
      }
    
    if (err < error) {
      converged <- TRUE
      break
    }
  }

on.exit()

if (progress)
  pb$tick(0, tokens = list(error = paste0(signif(err, 3), "/", signif(error, 3))))

if (verbose) {
  cat("Iterations performed:", i, "\n")
  cat("Converged:", converged, "\n")
}

if (!converged) {
  warning(
    "MDP solver did not converge after ",
    i,
    " iterations (delta = ",
    delta,
    ").",
    " Consider decreasing the 'discount' factor or increasing 'error' or 'N_max'."
  )
}

pi <- factor(pi, levels = seq_along(A), labels = A)

model$solution <- list(
  method = "value iteration (prioritized sweeping)",
  policy = list(data.frame(
    state = S,
    U = U,
    action = pi,
    row.names = NULL
  )),
  converged = converged,
  delta = delta,
  H = H,
  iterations = i
)
model
}


## Policy iteration

MDP_policy_iteration_inf_horizon <-
  function(model,
           N_max = 1000,
           k_backups = 10,
           U = NULL,
           matrix = TRUE,
           progress = TRUE,
           verbose = FALSE,
           ...) {
    if (progress) {
      pb <- my_progress_spinner(
        name = "solve_MDP",
        format_extra = paste0(
          " | changed actions: :changed_actions/",
          length(model$states),
          " | press esc/CTRL-C to terminate early"
        )
      )
      pb$tick(0, token = list(changed_actions = "-"))
    }
    
    S <- model$states
    
    if (is.null(U))
      U <- 0
    
    pi <- greedy_policy(q_values(model, U))
    U <- pi$U
    pi <- pi$action
    names(U) <- S
    names(pi) <- S
    
    if (progress)
      pb$tick(0, token = list(changed_actions = "-"))
    
    # return unconverged result when interrupted
    on.exit({
      warning("MDP solver did not converge (manual interrupt).")
      
      if (verbose) {
        cat("\nTerminated during iteration:", i, "\n")
      }
      
      model$solution <- list(
        method = "policy iteration",
        policy = list(data.frame(
          state = S,
          U = U,
          action = pi,
          row.names = NULL
        )),
        converged = converged,
        iterations = i
      )
      
      return(model)
    })
    
    changed_actions <- length(model$states)
    converged <- FALSE
    for (i in seq_len(N_max)) {
      if (progress)
        pb$tick(token = list(changed_actions = changed_actions))
      
      # evaluate to get U from pi
      policy_evaluation(
        model,
        pi,
        U,
        k_backups = k_backups,
        matrix = matrix,
        progress = FALSE
      )
      
      #if (progress)
      #  pb$tick(0, token = list(changed_actions = changed_actions))
      
      U_t_plus_1 <- bellman_update(model, U, return_Q = TRUE)
      pi_prime <- U_t_plus_1$pi
      U <- U_t_plus_1$U
      
      # account for random tie breaking using Q
      # if (all(pi == pi_prime)) {
      Q <- U_t_plus_1$Q
      changed_actions <- sum(Q[cbind(seq_len(nrow(Q)), pi)] != apply(Q, MARGIN = 1, max))
      
      if (verbose)
        cat(i, "| # of updated actions: ", changed_actions, "\n")
      
      
      if (changed_actions == 0) {
        converged <- TRUE
        break
      }
      
      pi <- pi_prime
    }
    
    #if (progress)
    #  pb$tick(0, token = list(changed_actions = changed_actions))
    
    # deregister on.exit
    on.exit()
    
    if (!converged) {
      warning(
        "MDP solver did not converge after ",
        i,
        " iterations.\n",
        " Consider decreasing the 'discount' factor or increasing 'N_max'."
      )
    }
    
    model$solution <- list(
      method = "policy iteration",
      policy = list(data.frame(
        state = S,
        U = U,
        action = pi,
        row.names = NULL
      )),
      converged = converged,
      iterations = i
    )
    model
  }
