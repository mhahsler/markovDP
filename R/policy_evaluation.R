#' Policy Evaluation
#'
#' Estimate the value function for a policy.
#' 
#' # Policy Evaluation using Bellman Operator
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
#' @param method used method: `"bellman"`, `"LP"`, or `"MC"`.
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
#' # compare the approx. value functions for the policies (we restrict
#' #    the number of backups for the random policy since it may not converge)
#' rbind(
#'   random = policy_evaluation(Maze, pi_random, k_backups = 100),
#'   manual = policy_evaluation(Maze, pi_manual),
#'   greedy = policy_evaluation(Maze, pi_greedy),
#'   optimal = policy_evaluation(Maze, pi_opt)
#' )
#' 
#' # use fist-visit Monte Carlo prediction with 100 episodes 
#' #   and a max horizon of 100
#' policy_evaluation(Maze, pi_opt, method = "MC", n = 100, horizon = 100)
#' @export
policy_evaluation <-
  function(model,
           pi = NULL,
           method = "bellman",
           ...,
           progress = TRUE,
           verbose = FALSE) {
  
    method <- match.arg(method, c("bellman", "LP", "MC")) 
    
    func <- get(paste0("policy_evaluation_", method))
    func(model, pi, ..., progress = progress, verbose = verbose)
}
