#' Policy Evaluation
#'
#' Estimate the value function for a policy applied to a 
#' model by repeatedly applying the Bellman operator.
#' 
#' The value function for a policy can be estimated (called policy evaluation)
#' by repeatedly applying the Bellman operator 
#' \deqn{v \leftarrow B_\pi(v)}
#' till convergence.
#' 
#' In each iteration, all state values are updated. In this implementation updating is
#' stopped when the largest state Bellman error is below a threshold.
#'
#' \deqn{||v_{k+1} - v_k||_\infty < \theta.}
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
#' maze_solved <- solve_MDP(Maze, method = "DP:VI")
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
#' pi_random <- random_policy(Maze, prob = c(up = .7, right = .1, down = .1, left = 0.1))
#' pi_random
#'
#' # 4. an improved policy based on one policy evaluation and
#' #   policy improvement step.
#' V <- policy_evaluation(Maze, pi_random)
#' Q <- Q_values(Maze, V)
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
    gamma <- model$discount %||% 1
    
    # start with all 0s if no previous U is given
    if (is.null(V)) {
      V <- V_zero(model)
    }
    
    # we cannot count more than integer.max
    if (k_backups > .Machine$integer.max) {
      k_backups <- .Machine$integer.max
      warning("Using the maximum number of backups of", k_backups)
    }

    # use r_pi and p_pi to make it faster and not do a complete Bellman update
    p_pi <- induced_transition_matrix(model, pi)
    r_pi <- induced_reward_matrix(model, pi)
    r_pi <- rowSums(p_pi * r_pi)
    
    for (i in seq_len(k_backups)) {
      if (progress)
        pb$tick()
      
      v <- V

      #V <- bellman_operator(model, pi, V)
      # do this faster directly
      V <- r_pi + gamma * p_pi %*% V
       
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
    
    V <- drop(V)
    names(V) <- S(model)
    V
  }




