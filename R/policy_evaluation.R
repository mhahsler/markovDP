#' Policy Evaluation
#'
#' Estimate the value function for a policy applied to a 
#' model by repeatedly applying the Bellman operator.
#' 
#' The Bellman operator updates a value function given the model defining
#' \eqn{T}, \eqn{\gamma} and \eqn{R}, and a policy
#' \eqn{\pi} by applying the Bellman equation as an update rule for each state:
#'
#' \deqn{U_{k+1}(s) =\sum_a \pi_{a|s} \sum_{s'} T(s' | s,a) [R(s,a) + \gamma U_k(s')]}
#'
#' A policy can be evaluated by applying the Bellman operator till convergence.
#' In each iteration, all states are updated. In this implementation updating is
#' stopped after`k_backups` iterations or after the
#' largest update
#'
#' \deqn{||U_{k+1} - U_k||_\infty < \theta.}
#'
#' @family MDP
#' @family policy
#' @author Michael Hahsler
#'
#' @param model an MDP problem specification.
#' @param U a vector with value function representing the state utilities
#'    (expected sum of discounted rewards from that point on).
#'    If `model` is a solved model, then the state
#'    utilities are taken from the solution.
#' @param pi a policy as a data.frame with at least columns for states and action.
#' @param k_backups number of look ahead steps used for approximate policy evaluation
#'    used by the policy iteration method. Set k_backups to `Inf` to only use
#'    \eqn{\theta} as the stopping criterion.
#' @param theta stop when the largest change in a state value is less
#'    than \eqn{\theta}.
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
#' u <- policy_evaluation(Maze, pi_random)
#' q <- q_values(Maze, U = u)
#' pi_greedy <- greedy_policy(q)
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
           pi,
           U = NULL,
           k_backups = 1000,
           theta = 1e-3,
           matrix = TRUE,
           progress = TRUE,
           verbose = FALSE) {
    if (!inherits(model, "MDP")) {
      stop("'model' needs to be of class 'MDP'.")
    }

    if (progress) {
      pb <- my_progress_spinner(name = "policy_evaluation")
      pb$tick(0)
    }
    
    S <- model$states
    A <- model$actions
    
    # note for the many iterations, it is better to get the complete matrices
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
        cat("done\n")
    }
    
    GAMMA <- model$discount

    if (is.data.frame(pi)) {
      pi <- pi$action
    }
    names(pi) <- S

    # start with all 0s if no previous U is given
    if (is.null(U)) {
      U <- rep(0, times = length(S))
    }
    names(U) <- S
    
    # we cannot count more than integer.max
    if (k_backups > .Machine$integer.max) {
      k_backups <- .Machine$integer.max
      warning("Using the maximum number of backups of", k_backups)
    }

    for (i in seq_len(k_backups)) {
      if (progress)
        pb$tick()
      
      v <- U
      
      # apply the Bellman operator
      if (matrix)
        U <- .QV_vec(S, pi, P, R, GAMMA, U)
      else
        U <- .QV_func_vec(S, pi, model, GAMMA, U)
      
      delta <- max(abs(v - U), na.rm = TRUE)

      if (verbose) {
        cat("Backup step", i, ": delta =", delta, "\n")
      }

      if (delta < theta) {
        break
      }
    }

    if (progress)
      pb$terminate()
    
    U
  }

### this is used if we already have a matrix
.policy_evaluation_int <-
  function(S,
           A,
           P,
           R,
           pi,
           GAMMA = 1,
           U = NULL,
           k_backups = 1000,
           theta = 1e-3,
           verbose = FALSE 
           ) {
   
    
    if (is.data.frame(pi)) {
      pi <- pi$action
    }
    names(pi) <- S
    
    # start with all 0s if no previous U is given
    if (is.null(U)) {
      U <- rep(0, times = length(S))
    }
    names(U) <- S
    
    # we cannot count more than integer.max
    if (k_backups > .Machine$integer.max) {
      k_backups <- .Machine$integer.max
      warning("Using the maximum number of backups of", k_backups)
    }
    
    for (i in seq_len(k_backups)) {
      v <- U
      
      # apply the Bellman operator
      U <- .QV_vec(S, pi, P, R, GAMMA, U)
      delta <- max(abs(v - U), na.rm = TRUE)
      
      if (verbose) {
        cat("Backup step", i, ": delta =", delta, "\n")
      }
      
      if (delta < theta) {
        break
      }
    }
    
    U
  }


#' @rdname policy_evaluation
#' @export
bellman_operator <- function(model, pi, U)
  .QV_func_vec(model$states, pi$action, model, model$discount %||% 1, U)


