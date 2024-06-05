# Solve MDPs using Dynamic Programming (Value and policy iteration)

# Calculate Q-Function from U
# the (optimal) state-action value function Q_s(a,k) is the expected total reward
# from stage k onward, if we choose a_k = a and then proceed optimally (given by U).
.QV <-
  function(s, a, P, R, GAMMA, U) {
    sum(P[[a]][s, ] * (R[[a]][[s]] + GAMMA * U), na.rm = TRUE)
  }
.QV_vec <- Vectorize(.QV, vectorize.args = c("s", "a"))


#' @rdname solve_MDP
#' @param U a vector with initial utilities used for each state. If
#'   `NULL`, then the default of a vector of all 0s is used.
#' @param N_max maximum number of iterations allowed to converge. If the
#'   maximum is reached then the non-converged solution is returned with a
#'   warning.
#' @param error value iteration: maximum error allowed in the utility of any state
#'   (i.e., the maximum policy loss) used as the termination criterion.
#' @param k_backups policy iteration: number of look ahead steps used for approximate policy evaluation
#'   used by the policy iteration method.
#' @export
solve_MDP_DP <- function(model,
                         method = "value_iteration",
                         horizon = NULL,
                         discount = NULL,
                         N_max = 1000,
                         error = 0.01,
                         k_backups = 10,
                         U = NULL,
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

  ret <- switch(method,
    value_iteration = {
      if (is.infinite(model$horizon)) {
        MDP_value_iteration_inf_horizon(model,
          error,
          N_max,
          U = U,
          verbose = verbose
        )
      } else {
        MDP_value_iteration_finite_horizon(model,
          horizon = model$horizon,
          U = U,
          verbose = verbose
        )
      }
    },
    policy_iteration = {
      if (is.infinite(model$horizon)) {
        MDP_policy_iteration_inf_horizon(model,
          N_max,
          k_backups,
          U = U,
          verbose = verbose
        )
      } else {
        stop("Method not implemented yet for finite horizon problems.")
      }
    },
    prioritized_sweeping = {
      if (is.infinite(model$horizon)) {
        MDP_PS_inf_horizon(model,
          error = error, N_max = N_max,
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
           verbose = FALSE) {
    if (!inherits(model, "MDP")) {
      stop("'model' needs to be of class 'MDP'.")
    }

    S <- model$states
    A <- model$actions
    P <- transition_matrix(model, sparse = TRUE)
    R <-
      reward_matrix(model, sparse = FALSE) ## Note sparse leaves it as a data.frame
    GAMMA <- model$discount

    horizon <- as.integer(horizon)

    if (is.null(U)) {
      U <- rep(0, times = length(S))
    }
    names(U) <- S

    policy <- vector(mode = "list", length = horizon)

    for (t in seq(horizon, 1)) {
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
           verbose = FALSE) {
    S <- model$states
    A <- model$actions
    P <- transition_matrix(model)
    R <- reward_matrix(model)
    GAMMA <- model$discount
    if (GAMMA < 1) {
      convergence_factor <- (1 - GAMMA) / GAMMA
    } else {
      convergence_factor <- 1
    }

    if (is.null(U)) {
      U <- rep(0, times = length(S))
    }
    names(U) <- S

    converged <- FALSE
    for (i in seq_len(N_max)) {
      if (verbose) {
        cat("Iteration:", i)
      }

      Qs <- outer(S, A, .QV_vec, P, R, GAMMA, U)
      m <- apply(Qs, MARGIN = 1, which.max.random)

      pi <- factor(m, levels = seq_along(A), labels = A)
      U_t_minus_1 <- Qs[cbind(seq_along(S), m)]

      delta <- max(abs(U_t_minus_1 - U))
      U <- U_t_minus_1

      if (verbose) {
        cat(" -> delta:", delta, "\n")
      }

      if (delta <= error * convergence_factor) {
        converged <- TRUE
        break
      }
    }

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
           verbose = FALSE) {
    S <- model$states
    A <- model$actions
    P <- transition_matrix(model)
    R <- reward_matrix(model)
    GAMMA <- model$discount
    if (GAMMA < 1) {
      convergence_factor <- (1 - GAMMA) / GAMMA
    } else {
      convergence_factor <- 1
    }

    if (is.null(U)) {
      U <- rep(0, times = length(S))
    }
    names(U) <- S

    pi <- rep(NA_integer_, times = length(S))

    ## State priorities. Set all priotities to error so we randomly pick one first.
    # chosen at least once.

    H <- rep(error, times = length(S))
    converged <- FALSE
    for (i in seq_len(N_max)) {
      if (verbose) {
        cat("Iteration", i)
      }

      # only update state with highest priority
      s <- which.max.random(H)

      Qs <- .QV_vec(s, A, P, R, GAMMA, U)
      m <- which.max.random(Qs)
      pi[s] <- m

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

      if (verbose) {
        cat(": state", s, " -> delta:", delta, "sum(H):", sum(H), "\n")
      }

      ### FIXME
      if (sum(H) < error) {
        converged <- TRUE
        break
      }
    }

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
      method = "value iteration",
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
           verbose = FALSE) {
    S <- model$states
    A <- model$actions
    P <- transition_matrix(model, sparse = TRUE)
    R <- reward_matrix(model, sparse = FALSE)
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

    converged <- FALSE
    for (i in seq_len(N_max)) {
      if (verbose) {
        cat("Iteration:", i, "\n")
      }

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
      if (all(Qs[cbind(seq_len(nrow(Qs)), pi)] == apply(Qs, MARGIN = 1, max))) {
        converged <- TRUE
        break
      }

      pi <- pi_prime
    }

    if (!converged) {
      warning(
        "MDP solver did not converge after ",
        i,
        " iterations.",
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
