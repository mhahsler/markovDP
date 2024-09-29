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

  # Note: This is used by gridworld_animate()
  dots <- list(...)
  if (!is.null(dots$n)) 
    iter_max <- dots$n
  
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
    message("No discount rate specified. Using .9 for the infinite horizon problem!")
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
          matrix = matrix,
          progress = progress,
          verbose = verbose,
          ...
        )
      } else {
        MDP_value_iteration_finite_horizon(model,
          horizon = model$horizon,
          U = U,
          matrix = matrix,
          progress = progress,
          verbose = verbose,
          ...
        )
      }
    },
    policy_iteration = {
      if (is.infinite(model$horizon)) {
        MDP_policy_iteration_inf_horizon(model,
          iter_max,
          k_backups,
          U = U,
          matrix = matrix,
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
        MDP_PS_inf_horizon(model,
          error = error, 
          iter_max = iter_max,
          U = U,
          matrix = matrix,
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
           matrix = TRUE,
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
    A <- model$actions
    
    if (matrix) {
      if (verbose)
        cat("Extracting matrices for R and T ...")
      R <- reward_matrix(model, sparse = NULL)
      if (progress)
        pb$tick(0)
      
      P <- transition_matrix(model, sparse = NULL)
      if (progress)
        pb$tick(0)
      if (verbose)
        cat("done.\n")
    }
    
    GAMMA <- model$discount
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
      
      if (matrix)
        Qs <- outer(S, A, .QV_vec, P, R, GAMMA, U)
      else
        Qs <- outer(S, A, .QV_func_vec, model, GAMMA, U)
      
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
      
    if (progress)
      pb$terminate()

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
           matrix = TRUE,
           progress = TRUE,
           verbose = FALSE,
           ...) {
    
    if (progress) {
      pb <- my_progress_spinner(name = "solve_MDP", 
                                format_extra = " | max delta U: :delta | press esc/CTRL-C to terminate early")
      pb$tick(0, token = list(delta = "-"))
    }
    
    S <- model$states
    A <- model$actions
    
    if (matrix) {
      if (verbose)
        cat("Extracting matrices for R and T ...")
      R <- reward_matrix(model, sparse = NULL)
      if (progress)
        pb$tick(0, token = list(delta = "-"))
      
      P <- transition_matrix(model, sparse = NULL)
      if (progress)
        pb$tick(0, token = list(delta = "-"))
      if (verbose)
        cat("done.\n")
    }
    
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
      if (matrix)
        Qs <- outer(S, A, .QV_vec, P, R, GAMMA, U)
      else
        Qs <- outer(S, A, .QV_func_vec, model, GAMMA, U)
      
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
##
## We initialize H(s) using the states reward so we start with the biggest 
## reward states and propagate the reward backwards.
MDP_PS_inf_horizon <-
  function(model,
           error,
           iter_max = 10000,
           U = NULL,
           init_H = "value_iteration",
           matrix = TRUE,
           progress = TRUE,
           verbose = FALSE,
           ...) {
    if (verbose)
      progress <- FALSE
    
    if (progress) {
      pb <- my_progress_spinner(name = "solve_MDP", ticks = "state updates",
                                format_extra = " | max priority: :error | press esc/CTRL-C to terminate early")
      pb$tick(0, token = list(error = "-"))
    }
    
    S <- model$states
    A <- model$actions
    
    if (matrix) {
      if (verbose)
        cat("Extracting matrices for R and T ...")
       
      R <- reward_matrix(model, sparse = NULL)
      if (progress) 
        pb$tick(0, token = list(error = "-"))
      
      P <- transition_matrix(model, sparse = NULL)
      if (progress) 
        pb$tick(0, token = list(error = "-"))
      
      if (verbose)
        cat("done\n")
    }
    
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
    
    # initialize with the maximum reward. This is a single value iteration pass.
    init_H <- match.arg(init_H, c("value_iteration", "random"))
     
    if (is.null(H)) {
      if (init_H == "value_iteration") {
        if (verbose)
          cat("initializing priority H using the Belmann error.\n\n")
        U <- bellman_update(model, U = U)
        H <- abs(U) + error + 1e-6
      } else { ### random
        if (verbose)
          cat("initializing priority H randomly greater than error.\n\n")
        H <- runif(length(model$states)) + error + 1e-6
      }
      
      # cat(round(H,2), "\n")
    }
    # may not find the optimal value if some priorities are 0 (See paper)

    names(U) <- S
    
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
        pb$tick(tokens = list(error = paste0(signif(err, 3), "/", 
                                             signif(error, 3))))

      # update state with highest priority (error) next
      s <- which.max.random(H)

      if (verbose)
        cat(i, "- updating state", s, sQuote(S[s]),"- action", pi[s], "-> ")
      
      if (matrix)
        Qs <- .QV_vec(s, A, P, R, GAMMA, U)
      else
        Qs <- .QV_func_vec(s, A, model, GAMMA, U)
      
      m <- which.max.random(Qs)
      pi[s] <- m
      delta <- abs(Qs[m] - U[s])
      
      if (verbose)
        cat(pi[s], "- U ", signif(U[s], 3), "->", Qs[m],"\n")
      
      U[s] <- Qs[m]
      
      # check for issues in pi
      if(any(is.na(pi)) || any(pi > length(A)) || any(pi < 1))
        stop("Problem with policy vector : ", pi)
      
      
      # This is Moore and Atkeson (1993) called PS in Li and Littman (2008)
      H[s] <- 0 # it will be updated below if it is a an ancestor

      # find states that can get us to s
      ancestor_prob <- transition_matrix(model, end.state = s, 
                                     simplify = TRUE, 
                                     sparse = TRUE)
      ancestor_ids <- Matrix::which(Matrix::rowSums(ancestor_prob, 
                                                    sparseResult = TRUE) > 0)
      
      # note this also increases the priority for the state we just updated from
      for (ancestor_id in ancestor_ids) 
        H[ancestor_id] <- max(H[ancestor_id], delta * max(ancestor_prob[ancestor_id, ]))
      
      err <- max(H)
      
      if (verbose) {
        cat("    updating H for states:", paste(ancestor_ids, sQuote(S[ancestor_ids]) , collapse = ", "),"\n")
      }
        
      
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
      pb <- my_progress_spinner(name = "solve_MDP", 
                                format_extra = paste0(" | changed actions: :changed_actions/",
                                                      length(model$states),
                                                      " | press esc/CTRL-C to terminate early"))
      pb$tick(0, token = list(changed_actions = "-"))
    }
      
    S <- model$states
    A <- model$actions
    
    if (matrix) {
      if (verbose)
        cat("Extracting matrices for R and T ...")
      R <- reward_matrix(model, sparse = NULL)
      if (progress)
        pb$tick(0, token = list(changed_actions = "-"))
      
      P <- transition_matrix(model, sparse = NULL)
      if (progress)
        pb$tick(0, token = list(changed_actions = "-"))
      if (verbose)
        cat("done.\n")
    }
    
    GAMMA <- model$discount

    if (is.null(U)) {
      U <- rep(0, times = length(S))
      pi <- random_policy(model)$action
    } else {
      # get greedy policy for a given U
      if (matrix)
        Qs <- outer(S, A, .QV_vec, P, R, GAMMA, U)
      else
        Qs <- outer(S, A, .QV_func_vec, model, GAMMA, U)
      
      m <- apply(Qs, MARGIN = 1, which.max.random)
      pi <- factor(m, levels = seq_along(A), labels = A)
    }

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
      if (matrix)
        U <- .policy_evaluation_int(S, A, P, R, pi, GAMMA = GAMMA, U =  U, k_backups = k_backups)
      else
        U <-
          policy_evaluation(model, pi, U, k_backups = k_backups, matrix = matrix, progress = FALSE)
        
      if (progress)
        pb$tick(0, token = list(changed_actions = changed_actions))
      
      # get greedy policy for U
      if (matrix)
        Qs <- outer(S, A, .QV_vec, P, R, GAMMA, U)
      else
        Qs <- outer(S, A, .QV_func_vec, model, GAMMA, U)

      if (progress)
        pb$tick(0, token = list(changed_actions = changed_actions))
      
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
