# TODO: deal with available actions for states actions(s)

#' Solve an MDP Problem
#'
#' Implementation of value iteration, modified policy iteration and other
#' methods based on reinforcement learning techniques to solve finite
#' state space MDPs.
#'
#'
#' Several solvers are available.
#'
#' ## Dynamic Programming
#' Implemented are the following dynamic programming methods (following
#' Russell and Norvig, 2010):
#'
#' * **Modified Policy Iteration**
#' starts with a random policy and iteratively performs
#' a sequence of
#'   1. approximate policy evaluation (estimate the value function for the
#' current policy using `k_backups` and function [`policy_evaluation()`], and
#'   2. policy improvement (calculate a greedy policy given the value function).
#' The algorithm stops when it converges to a stable policy (i.e., no changes
#' between two iterations).
#'
#' * **Value Iteration** starts with
#'   an arbitrary value function (by default all 0s) and iteratively
#'   updates the value function for each state using the Bellman equation.
#'   The iterations
#'   are terminated either after `N_max` iterations or when the solution converges.
#'   Approximate convergence is achieved
#'   for discounted problems (with \eqn{\gamma < 1})
#'   when the maximal value function change for any state \eqn{\delta} is
#'   \eqn{\delta \le error (1-\gamma) / \gamma}. It can be shown that this means
#'   that no state value is more than
#'   \eqn{error} from the value in the optimal value function. For undiscounted
#'   problems, we use \eqn{\delta \le error}.
#'
#'   The greedy policy
#'   is calculated from the final value function. Value iteration can be seen as
#'   policy iteration with truncated policy evaluation.
#'
#' Note that the policy converges earlier than the value function.
#'
#' ## Linear Programming
#' The following linear programming formulation is implemented. For the
#' optimal value function, the Bellman equation holds:
#'
#' \deqn{
#'  V^*(s) = \max_{a \in A}\sum_{s' \in S} T(s, a, s') [ R(s, a, s') + \gamma V^*(s')]\; \forall a\in A, s \in S
#' }
#'
#' We can find the optimal value function by solving the following linear program:
#' \deqn{\text{min} \sum_{s\in S} V(s)}
#' subject to
#' \deqn{V(s) \ge \sum_{s' \in S} T(s, a, s') [R(s, a, s') + \gamma V(s')],\; \forall a\in A, s \in S
#' }
#'
#' Note:
#' * The discounting factor has to be strictly less than 1.
#' * Additional parameters to to `solve_MDP` are passed on to [lpSolve::lp()].
#' * We use the solver in
#'   [lpSolve::lp()] which requires all decision variables (state values) to be non-negative.
#'   To ensure this, for negative rewards, all rewards as shifted so the
#'   smallest reward is
#'   0. This does not change the optimal policy.
#'
#' ## Temporal Difference Control
#'
#' Implemented are the following temporal difference control methods
#' described in Sutton and Barto (2020).
#' Note that the MDP transition and reward models are only used to simulate
#' the environment for these reinforcement learning methods.
#' The algorithms use a step size parameter \eqn{\alpha} (learning rate) for the
#' updates and the exploration parameter \eqn{\epsilon} for
#' the \eqn{\epsilon}-greedy policy.
#'
#' If the model has absorbing states to terminate episodes, then no maximal episode length
#' (`horizon`) needs to
#' be specified. To make sure that the algorithm does finish in a reasonable amount of time,
#' episodes are stopped after 10,000 actions with a warning. For models without absorbing states,
#' a episode length has to be specified via `horizon`.
#'
#' * **Q-Learning** is an off-policy temporal difference method that uses
#'    an \eqn{\epsilon}-greedy behavior policy and learns a greedy target
#'    policy.
#'
#' * **Sarsa** is an on-policy method that follows and learns
#'    an \eqn{\epsilon}-greedy policy. The final \eqn{\epsilon}-greedy policy
#'    is converted into a greedy policy.
#'
#' * **Expected Sarsa**: We implement an on-policy version that uses
#'   the expected value under the current policy for the update.
#'   It moves deterministically in the same direction as Sarsa
#'   moves in expectation. Because it uses the expectation, we can
#'   set the step size \eqn{\alpha} to large values and even 1.
#'
#' @family solver
#' @family MDP
#'
#' @param model an MDP problem specification.
#' @param method string; one of the following solution methods: `'value_iteration'`,
#'   `'policy_iteration'`, `'lp'`, `'q_learning'`, `'sarsa'`, or `'expected_sarsa'`.
#' @param horizon an integer with the number of epochs for problems with a
#'   finite planning horizon. If set to `Inf`, the algorithm continues
#'   running iterations till it converges to the infinite horizon solution. If
#'   `NULL`, then the horizon specified in `model` will be used.
#' @param discount discount factor in range \eqn{(0, 1]}. If `NULL`, then the
#'   discount factor specified in `model` will be used.
#' @param verbose logical, if set to `TRUE`, the function provides the
#'   output of the solver in the R console.
#' @param ... further parameters are passed on to the solver function.
#'
#' @return `solve_MDP()` returns an object of class POMDP which is a list with the
#'   model specifications (`model`), the solution (`solution`).
#'   The solution is a list with the elements:
#'   - `policy` a list representing the policy graph. The list only has one element for converged solutions.
#'   - `converged` did the algorithm converge (`NA`) for finite-horizon problems.
#'   - `delta` final \eqn{\delta} (value iteration and infinite-horizon only)
#'   - `iterations` number of iterations to convergence (infinite-horizon only)
#'
#' @author Michael Hahsler
#' @references
#' Russell, S., Norvig, P. (2021). Artificial Intelligence: A Modern Approach.
#' Fourth edition. Prentice Hall.
#'
#' Sutton, R. S., Barto, A. G. (2020). Reinforcement Learning: An Introduction.
#' Second edition. The MIT Press.
#'
#' @examples
#' data(Maze)
#' Maze
#'
#' # use value iteration
#' maze_solved <- solve_MDP(Maze, method = "value_iteration")
#' maze_solved
#' policy(maze_solved)
#'
#' # plot the value function U
#' plot_value_function(maze_solved)
#'
#' # Gridworld solutions can be visualized
#' gridworld_plot(maze_solved)
#'
#' # Use linear programming
#' maze_solved <- solve_MDP(Maze, method = "lp")
#' maze_solved
#' policy(maze_solved)
#'
#' # use modified policy iteration
#' maze_solved <- solve_MDP(Maze, method = "policy_iteration")
#' policy(maze_solved)
#'
#' # finite horizon
#' maze_solved <- solve_MDP(Maze, method = "value_iteration", horizon = 3)
#' policy(maze_solved)
#' gridworld_plot(maze_solved, epoch = 1)
#' gridworld_plot(maze_solved, epoch = 2)
#' gridworld_plot(maze_solved, epoch = 3)
#'
#' # create a random policy where action n is very likely and approximate
#' #  the value function. We change the discount factor to .9 for this.
#' Maze_discounted <- Maze
#' Maze_discounted$discount <- .9
#' pi <- random_policy(Maze_discounted,
#'   prob = c(n = .7, e = .1, s = .1, w = 0.1)
#' )
#' pi
#'
#' # compare the utility function for the random policy with the function for the optimal
#' #  policy found by the solver.
#' maze_solved <- solve_MDP(Maze)
#'
#' policy_evaluation(Maze, pi, k_backup = 100)
#' policy_evaluation(Maze, policy(maze_solved), k_backup = 100)
#'
#' # Note that the solver already calculates the utility function and returns it with the policy
#' policy(maze_solved)
#'
#' # Learn a Policy using Q-Learning
#' maze_learned <- solve_MDP(Maze, method = "q_learning", N = 100)
#' maze_learned
#'
#' maze_learned$solution
#' policy(maze_learned)
#' plot_value_function(maze_learned)
#' gridworld_plot(maze_learned)
#' @export
solve_MDP <- function(model, method = "value", ...) {
  methods_DP <- c("value_iteration", "policy_iteration")
  methods_LP <- c("lp")
  methods_TD <- c("sarsa", "q_learning", "expected_sarsa")

  if (!inherits(model, "MDP")) {
    stop("x needs to be a MDP!")
  }

  method <- match.arg(method, c(methods_DP, methods_LP, methods_TD))

  if (method %in% methods_DP) {
    return(solve_MDP_DP(model, method, ...))
  }
  if (method %in% methods_LP) {
    return(solve_MDP_LP(model, method, ...))
  }
  if (method %in% methods_TD) {
    return(solve_MDP_TD(model, method, ...))
  }

  # we should not get here!
  stop("Unknown method!")
}
