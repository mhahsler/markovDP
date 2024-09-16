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
#' @export
solve_MDP_DP <- function(model,
                         method = "value_iteration",
                         horizon = NULL,
                         discount = NULL,
                         iter_max = 1000,
                         error = 0.01,
                         k_backups = 10,
                         U = NULL,
                         continue = FALSE,
                         progress = TRUE,
                         verbose = FALSE) {
  if (!inherits(model, "MDP")) {
    stop("'model' needs to be of class 'MDP'.")
  }

  methods <- c("value_iteration", "policy_iteration", "prioritized_sweeping")
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
    message("No discount rate specified. Using .9!")
    model$discount <- .9
  }

  if (continue) {
    if (is.null(model$solution$policy[[1]]$U))
      stop("model solution does not contain a U vector to continue from!")
    U <- model$solution$policy[[1]]$U
  }
  
  ret <- switch(method,
    value_iteration = {
      if (is.infinite(model$horizon)) {
        MDP_value_iteration_inf_horizon(model,
          error,
          iter_max,
          U = U,
          progress = progress,
          verbose = verbose
        )
      } else {
        MDP_value_iteration_finite_horizon(model,
          horizon = model$horizon,
          U = U,
          progress = progress,
          verbose = verbose
        )
      }
    },
    policy_iteration = {
      if (is.infinite(model$horizon)) {
        MDP_policy_iteration_inf_horizon(model,
          iter_max,
          k_backups,
          U = U,
          progress = progress,
          verbose = verbose
        )
      } else {
        stop("Method not implemented yet for finite horizon problems.")
      }
    },
    prioritized_sweeping = {
      if (is.infinite(model$horizon)) {
        MDP_PS_inf_horizon(model,
          error = error, 
          N_max = iter_max,
          U = U,
          progress = progress,
          verbose = verbose
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
           verbose = FALSE) {
    if (!inherits(model, "MDP")) {
      stop("'model' needs to be of class 'MDP'.")
    }

    S <- model$states
    A <- model$actions
    P <- transition_matrix(model, sparse = TRUE)
    R <- reward_matrix(model, sparse = TRUE)
    GAMMA <- model$discount

    horizon <- as.integer(horizon)

    if (is.null(U)) {
      U <- rep(0, times = length(S))
    }
    names(U) <- S

    policy <- vector(mode = "list", length = horizon)

    if (progress) {
      pb <- my_progress_bar(horizon, name = "solve_MDP") 
      pb$tick(0)
    }
      
    for (t in seq(horizon, 1)) {
      if (progress)
        pb$tick()
      
      
      if (verbose) {
        cat("Iteration for t = ", t)
      }

      Qs <- outer(S, A, .QV_vec, P, R, GAMMA, U)
      m <- apply(Qs, MARGIN = 1, which.max.random)

      pi <- factor(m, levels = seq_along(A), labels = A)
      U <- Qs[cbind(seq_along(S), m)]

      policy[[t]] <- data.frame(
        state = S,
        U = U,
        action = pi,
        row.names = NULL
      )
    }

    model$solution <- list(
      policy = policy,
      converged = NA
    )
    model
  }

MDP_value_iteration_inf_horizon <-
  function(model,
           error,
           N_max = 1000,
           U = NULL,
           progress = TRUE,
           verbose = FALSE) {
    S <- model$states
    A <- model$actions
    P <- transition_matrix(model, sparse = TRUE)
    R <- reward_matrix(model, sparse = TRUE)
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

    if (progress)
      pb <- my_progress_spinner(name = "solve_MDP", 
                                format_extra = " | max delta U: :delta | press esc/CTRL-C to terminate early")
    
    # return unconverged result when interrupted
    on.exit({ 
      warning("MDP solver did not converge (manual interrupt).")
      
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
                  list(delta = paste0(signif(delta, 3), "/", 
                                       signif(convergence_limit, 3))))
      
      if (verbose && !progress) {
        cat("Iteration:", i)
      }

      # Find the best action for each state given U
      Qs <- outer(S, A, .QV_vec, P, R, GAMMA, U)
      m <- apply(Qs, MARGIN = 1, which.max.random)
      pi <- factor(m, levels = seq_along(A), labels = A)
      U_t_plus_1 <- Qs[cbind(seq_along(S), m)]

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
                list(delta = paste0(signif(delta, 3), "/", 
                                    signif(convergence_limit, 3))))

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
MDP_PS_inf_horizon <-
  function(model,
           error,
           N_max = 1000,
           U = NULL,
           progress = TRUE,
           verbose = FALSE) {
    S <- model$states
    A <- model$actions
    P <- transition_matrix(model, sparse = TRUE)
    R <- reward_matrix(model, sparse = TRUE)
    GAMMA <- model$discount

    if (is.null(U)) {
      U <- rep(0, times = length(S))
      # random policy
      pi <- sample(seq_along(A), size = length(S), replace = TRUE)
      H <- rep(.reward_range(model)[2] *.1, times = length(S))
    } else {
      pi <- model$solution$policy$action
      if (is.null(pi))
        pi <- sample(seq_along(A), size = length(S), replace = TRUE)
      
      H <- model$solution$H
      if (is.null(H))
        H <- rep(.reward_range(model)[2] *.1, times = length(S))
      
      
    }
    names(U) <- S

    ## State priorities. Set all priorities to error so we randomly pick one first.
    # chosen at least once.

    if (progress)
      pb <- my_progress_spinner(name = "solve_MDP", ticks = "state updates",
                                format_extra = " | max priority: :error | press esc/CTRL-C to terminate early")
    
    # return unconverged result when interrupted
    on.exit({ 
      warning("MDP solver did not converge (manual interrupt).")
      
      if (verbose) {
        cat("\nTerminated during iteration:", i, "\n")
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
      return(model)
    })
    
    
    err <- sum(H)
    converged <- FALSE
    for (i in seq_len(N_max)) {
      if (progress)
        pb$tick(tokens = list(error = paste0(signif(err, 3), "/", 
                                             signif(error, 3))))

      # FIXME: Maybe we need a sweep first or it will just randomly try states
      #        till it finds one with a different reward (if step cost is 0).
      # update state with highest error
      s <- which.max.random(H)

      Qs <- .QV_vec(s, A, P, R, GAMMA, U)
      m <- which.max.random(Qs)
      pi[s] <- m

      # check for issues in pi
      if(any(is.na(pi)) || any(pi > 4) || any(pi < 1))
        stop("Problem with policy vector : ", pi)
      
      delta <- abs(Qs[m] - U[s])
      U[s] <- Qs[m]

      # update priority for all states
      for (ss in seq_along(S)) {
        if (ss == s) {
          H[ss] <- delta * max(sapply(A, FUN = function(a) P[[a]][s, ss]))
        } else {
          H[ss] <- max(H[ss], delta * max(sapply(A, FUN = function(a) P[[a]][s, ss])))
        }
      }


      # pi <- factor(m, levels = seq_along(A), labels = A)
      # U_t_minus_1 <- Qs[cbind(seq_along(S), m)]
      #
      # delta <- max(abs(U_t_minus_1 - U))
      # U <- U_t_minus_1

      err <- max(H)
      
      if (verbose && !progress) {
        cat(
          "Try ",
          format(i, width = 6) ,
          ": state",
          format(model$states[s], width = 10),
          "-> delta:",
          format(delta, width = 6),
          "max(H):",
          format(err, width = 6),
          "\n"
        )
      }
      
      ### FIXME
      if (err < error) {
        converged <- TRUE
        break
      }
    }

    on.exit()
    
    if (progress)
      pb$tick(0, tokens = list(error = paste0(signif(err, 3), "/", 
                                           signif(error, 3)))) 
    
    if (verbose) {
      cat("Iterations needed:", i, "\n")
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
           progress = TRUE,
           verbose = FALSE) {
    S <- model$states
    A <- model$actions
    P <- transition_matrix(model, sparse = TRUE)
    R <- reward_matrix(model, sparse = TRUE)
    GAMMA <- model$discount

    if (is.null(U)) {
      U <- rep(0, times = length(S))
      pi <- random_policy(model)$action
    } else {
      # get greedy policy for a given U
      Qs <- outer(S, A, .QV_vec, P, R, GAMMA, U)
      m <- apply(Qs, MARGIN = 1, which.max.random)
      pi <- factor(m, levels = seq_along(A), labels = A)
    }

    names(U) <- S
    names(pi) <- S
    
    if (progress)
      pb <- my_progress_spinner(name = "solve_MDP", format_extra = " | changed actions: :changed_actions | press esc/CTRL-C to terminate early")
    
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
      U <-
        policy_evaluation(model, pi, U, k_backups = k_backups)

      # get greedy policy for U
      Qs <- outer(S, A, .QV_vec, P, R, GAMMA, U)

      # non-randomizes
      m <- apply(Qs, MARGIN = 1, which.max.random)
      pi_prime <- factor(m, levels = seq_along(A), labels = A)

      # account for random tie breaking
      # if (all(pi == pi_prime)) {
      changed_actions <- sum(Qs[cbind(seq_len(nrow(Qs)), pi)] != apply(Qs, MARGIN = 1, max))
      if (changed_actions == 0) {
        converged <- TRUE
        break
      }

      pi <- pi_prime
    }

    if (progress)
      pb$tick(0, token = list(changed_actions = changed_actions))
    
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
