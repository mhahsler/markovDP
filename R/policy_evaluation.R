#' Policy Evaluation
#'
#' Estimate the value function for a policy applied to a 
#' model by repeatedly applying the Bellman operator.
#' 
#' The value function for a policy can be estimated (called policy evaluation)
#' by repeatedly applying the Bellman operator \eqn{B_\pi} till convergence.
#' In each iteration, all state values are updated. In this implementation updating is
#' stopped when the largest state Bellman error is below a threshold.
#'
#' \deqn{||V_{k+1} - V_k||_\infty < \theta.}
#'
#' Or if `k_backups` iterations have been completed.
#'
#' @family MDP
#' @family policy
#' @author Michael Hahsler
#'
#' @param model an MDP problem specification.
#' @param V a vector with estimated state values representing a value function.
#'    If `model` is a solved model, then the state
#'    values are taken from the solution.
#' @param pi a policy as a data.frame with at least columns for states and action. If `NULL`,
#'     then the policy in model is used.
#' @param k_backups number of look ahead steps used for approximate policy evaluation
#'    used by the policy iteration method. Set k_backups to `Inf` to only use
#'    \eqn{\theta} as the stopping criterion.
#' @param theta stop when the largest state Bellman error (\eqn{\delta = V_{k+1} - V}) 
#'    is less than \eqn{\theta}.
#' @param matrix logical; if `TRUE` then matrices for the transition model and 
#'    the reward function are taken from the model first. This can be slow if functions 
#'    need to be converted or not fit into memory if the models are large. If these
#'    components are already matrices, then this is very fast. For `FALSE`, the
#'    transition probabilities and the reward is extracted when needed. This is slower, 
#'    but removes the time to calculate the matrices and it saves memory.
#' @param progress logical; show a progress bar with estimated time for completion.
#' @param verbose logical; should progress and approximation errors be printed.
#' @return a vector with (approximate)
#'    state values (U).
#'
#' @references
#' Sutton, R. S., Barto, A. G. (2020). Reinforcement Learning: An Introduction.
#' Second edition. The MIT Press.
#'
#' @examples
#' data(Maze)
#' Maze
#'
#' # create several policies:
#' # 1. optimal policy using value iteration
#' maze_solved <- solve_MDP(Maze, method = "value_iteration")
#' pi_opt <- policy(maze_solved)
#' pi_opt
#'
#' # 2. a manual policy (go up and in some squares to the right)
#' acts <- rep("up", times = length(Maze$states))
#' names(acts) <- Maze$states
#' acts[c("s(1,1)", "s(1,2)", "s(1,3)")] <- "right"
#' pi_manual <- manual_policy(Maze, acts)
#' pi_manual
#'
#' # 3. a random policy
#' set.seed(1234)
#' pi_random <- random_policy(Maze)
#' pi_random
#'
#' # 4. an improved policy based on one policy evaluation and
#' #   policy improvement step.
#' V <- policy_evaluation(Maze, pi_random)
#' Q <- q_values(Maze, V)
#' pi_greedy <- greedy_policy(Q)
#' pi_greedy
#'
#' #' compare the approx. value functions for the policies (we restrict
#' #'    the number of backups for the random policy since it may not converge)
#' rbind(
#'   random = policy_evaluation(Maze, pi_random, k_backups = 100),
#'   manual = policy_evaluation(Maze, pi_manual),
#'   greedy = policy_evaluation(Maze, pi_greedy),
#'   optimal = policy_evaluation(Maze, pi_opt)
#' )
#' @export
policy_evaluation <-
  function(model,
           pi = NULL,
           V = NULL,
           k_backups = 1000L,
           theta = 1e-3,
           matrix = TRUE,
           progress = TRUE,
           verbose = FALSE) {
    if (!inherits(model, "MDP")) {
      stop("'model' needs to be of class 'MDP'.")
    }

    if (is.null(pi))
      pi <- policy(model)
    
    if (progress) {
      pb <- my_progress_spinner(name = "policy_evaluation")
      pb$tick(0)
    }
    
    S <- S(model)
    A <- A(model)
    
    # note for the many iterations, it is better to get the complete matrices
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
    
    if (is.data.frame(pi)) {
      pi <- pi$action
    }
    names(pi) <- S

    # start with all 0s if no previous U is given
    if (is.null(V)) {
      V <- rep(0, times = length(S))
    }
    names(V) <- S
    
    # we cannot count more than integer.max
    if (k_backups > .Machine$integer.max) {
      k_backups <- .Machine$integer.max
      warning("Using the maximum number of backups of", k_backups)
    }

    for (i in seq_len(k_backups)) {
      if (progress)
        pb$tick()
      
      v <- V
      V <- bellman_operator(model, pi, V)
      delta <- max(abs(v - V), na.rm = TRUE)

      if (verbose) {
        cat("Backup step", i, ": delta =", delta, "\n")
      }

      if (delta < theta) {
        break
      }
    }

    if (progress)
      pb$terminate()
    
    V
  }




