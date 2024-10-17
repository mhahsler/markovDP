#' Solve an MDP Problem
#'
#' Implementation of value iteration, modified policy iteration and other
#' methods based on reinforcement learning techniques to solve finite
#' state space MDPs.
#'
#'
#' Several solvers are available. Note that some solvers are only implemented
#' for finite-horizon problems.
#'
#' ## Dynamic Programming
#'
#' Implemented are the following dynamic programming methods (following
#' Russell and Norvig, 2010):
#'
#' * **Modified Policy Iteration** (Howard 1960; Puterman and Shin 1978)
#' starts with a random policy and iteratively performs
#' a sequence of
#'   1. approximate policy evaluation (estimate the value function for the
#' current policy using `k_backups` and function [`policy_evaluation()`], and
#'   2. policy improvement (calculate a greedy policy given the value function).
#' The algorithm stops when it converges to a stable policy (i.e., no changes
#' between two iterations).
#'
#' * **Value Iteration** (Bellman 1957) starts with
#'   an arbitrary value function (by default all 0s) and iteratively
#'   updates the value function for each state using the Bellman equation.
#'   The iterations
#'   are terminated either after `n` iterations or when the solution converges.
#'   Approximate convergence is achieved
#'   for discounted problems (with \eqn{\gamma < 1})
#'   when the maximal value function change for any state \eqn{\delta} is
#'   \eqn{\delta \le error (1-\gamma) / \gamma}. It can be shown that this means
#'   that no state value is more than
#'   \eqn{error} from the value in the optimal value function. For undiscounted
#'   problems, we use \eqn{\delta \le error}.
#'
#'   A greedy policy
#'   is calculated from the final value function. Value iteration can be seen as
#'   policy iteration with truncated policy evaluation.
#'
#' * **Prioritized Sweeping** (Moore and Atkeson, 1993; Andre et al., 1997; Li and Littman, 2008)
#'   approximate the optimal value
#'   function by iteratively adjusting one state at a time. The state to be updated is chosen
#'   depending on its priority which reflects how much a state value may change
#'   given the most recently updated other states that can be directly reached via an action.
#'   This update order often lead to faster convergence compared
#'   to sweeping the whole state state in regular value iteration.
#'
#'   We implement the two priority update strategies described as __PS__ and
#'   __GenPS__ by Li and Littman.
#'
#'   * __PS__ (Moore and Atkeson, 1993) updates the priority of a state \eqn{H(s)}
#'      using:
#'      \deqn{
#'        \forall{s \in S}: H_{t+1}(s)  \leftarrow \begin{cases}
#'          \max(H_{t}(s), \Delta_t \max_{a \in A}(T(s_t|s,a)) \text{ for } s \ne s_{t+1} \\
#'          \Delta_t \max_{a \in A}(T(s_t|s,a) \text{ for } s = s_{t+1}
#'          \end{cases}
#'      }
#'
#'      where \eqn{\Delta_t = |V_{t+1}(s_t) - V_t(s_t)| = |E(s_t; V_{t+1})|}, i.e.,
#'      the Bellman error for the updated state.
#'
#'   * __GenPS__ (Andre et al., 1997) updates all state priorities using their
#'      current Bellman error:
#'
#'      \deqn{\forall{s \in S}: H_{t+1}(s) \leftarrow |E(s; V_{t+1})|}
#'
#'      where \eqn{E(s; V_{t+1}) = \max_{a \in A} \left[R(s,a) + \gamma \sum_{s \in S} T(s'|s,a) V(s')\right] - V(s)}.
#'      is a state's Bellman error.
#'
#'   The update method can be chosen using the additional parameter `H_update`
#'   as the character string `"PS_random"`, `"PS_error"` or `"GenPS"`.
#'   The default is `H_update = "GenPS"`. For PS, random means that the
#'   priority vector is initialized with random values (larger than 0),
#'   and error means they are initialized with the Bellman error as in
#'   GenPS. However, this requires one complete sweep over all states.
#'
#'   This implementation stops updating when the largest priority values
#'   over all states is less than the specified `error`.
#'
#'   Since the algorithm does not sweep through the whole state space for each
#'   iteration, `n` is converted into an equivalent number of state updates
#'   \eqn{n = n |S|}.
#'
#' Note that policies converge earlier than value functions.
#'
#' ## Linear Programming
#'
#' The following linear programming formulation (Manne 1960) is implemented. For the
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
#'   0. This transformation does not change the optimal policy.
#'
#'
#' ## Monte Carlo Control
#'
#' Monte Carlo control simulates a whole episode using the current behavior
#' policy and then updates the target policy before simulating the next episode.
#' Implemented are the following temporal difference control methods
#' described in Sutton and Barto (2020).
#'
#' * **Monte Carlo Control with exploring Starts** uses the same greedy policy for
#' behavior and target (on-policy). To make sure all states/action pairs are
#' explored, it uses exploring starts meaning that new episodes are started at a randomly
#' chosen state using a randomly chooses action.
#'
#' * **On-policy Monte Carlo Control** uses for behavior and as the target policy
#' an epsilon-greedy policy.
#'
#' * **Off-policy Monte Carlo Control** uses for behavior an arbitrary policy
#' (we use an epsilon-greedy policy) and learns a greedy policy using
#' importance sampling.
#'
#' ## Temporal Difference Control
#'
#' Implemented are several temporal difference control methods
#' described in Sutton and Barto (2020).
#' Note that the MDP transition and reward models are used
#' for these reinforcement learning methods only to sample from
#' the environment.
#' The algorithms use a step size parameter \eqn{\alpha} (learning rate) for the
#' updates and the exploration parameter \eqn{\epsilon} for
#' the \eqn{\epsilon}-greedy behavior policy.
#'
#' If the model has absorbing states to terminate episodes, then no maximal episode length
#' (`horizon`) needs to
#' be specified. To make sure that the algorithm does finish in a reasonable amount of time,
#' episodes are stopped after 1000 actions (with a warning). For models without absorbing states,
#' the episode length has to be specified via `horizon`.
#'
#' * **Q-Learning** (Watkins and Dayan 1992) is an off-policy temporal difference method that uses
#'    an \eqn{\epsilon}-greedy behavior policy and learns a greedy target
#'    policy.
#'
#' * **Sarsa** (Rummery and Niranjan 1994) is an on-policy method that follows and learns
#'    an \eqn{\epsilon}-greedy policy. The final \eqn{\epsilon}-greedy policy
#'    is converted into a greedy policy.
#'
#' * **Expected Sarsa** (R. S. Sutton and Barto 2018). We implement an on-policy version that uses
#'   the expected value under the current policy for the update.
#'   It moves deterministically in the same direction as Sarsa would
#'   move in expectation. Because it uses the expectation, we can
#'   set the step size \eqn{\alpha} to large values and 1 is common.
#'
#' ## Planning by Sampling
#'
#' A simple, not very effective, planning method proposed by Sutton and Barto (2020) in Chapter 8.
#'
#' * **Random-sample one-step tabular Q-planning** randomly selects a
#' state/action pair and samples the resulting reward and next state from
#' the model. This
#' information is used to update a single Q-table value.
#'
#' @family solver
#' @family MDP
#'
#' @param model an MDP problem specification.
#' @param method string; one of the following solution methods:
#'    `'value_iteration'`,
#'    `'policy_iteration'`,
#'    `'lp'`,
#'    `'q_learning'`,
#'    `'sarsa'`,
#'    `'expected_sarsa'`,
#'    `'MC_exploring_starts'`,
#'    `'MC_on_policy'`,
#'    `'MC_off_policy',
#'    `'q_planning'`.
#' @param horizon an integer with the number of epochs for problems with a
#'   finite planning horizon. If set to `Inf`, the algorithm continues
#'   running iterations till it converges to the infinite horizon solution. If
#'   `NULL`, then the horizon specified in `model` will be used.
#' @param discount discount factor in range \eqn{(0, 1]}. If `NULL`, then the
#'   discount factor specified in `model` will be used.
#' @param progress logical; show a progress bar with estimated time for completion.
#' @param verbose logical, if set to `TRUE`, the function provides the
#'   output of the solver in the R console.
#' @param ... further parameters are passed on to the solver function.
#'
#' @return `solve_MDP()` returns an object of class MDP which is a list with the
#'   model specifications (`model`), the solution (`solution`).
#'   The solution is a list with the elements:
#'   - `policy` a list representing the policy graph. The list only has one
#'      element for converged solutions.
#'   - `converged` did the algorithm converge (`NA`) for finite-horizon problems.
#'   - `delta` final \eqn{\delta} (value iteration and infinite-horizon only)
#'   - `iterations` number of iterations to convergence (infinite-horizon only)
#'
#' @author Michael Hahsler
#' @references
#' Andre, D., Friedman, N., and Parr, R. 1997. "Generalized prioritized sweeping." In Advances in Neural Information Processing Systems 10, pp. 1001-1007. [NeurIPS Proceedings](https://proceedings.neurips.cc/paper_files/paper/1997/file/7b5b23f4aadf9513306bcd59afb6e4c9-Paper.pdf)
#'
#' Bellman, Richard. 1957. "A Markovian Decision Process." Indiana University Mathematics Journal 6: 679-84. [https://www.jstor.org/stable/24900506](https://www.jstor.org/stable/24900506).
#'
#' Howard, R. A. 1960. "Dynamic Programming and Markov Processes." Cambridge, MA: MIT Press.
#'
#' Li, Lihong, and Michael Littman. 2008. "Prioritized Sweeping Converges to the Optimal Value Function." DCS-TR-631. Rutgers University. \doi{10.7282/T3TX3JSX}
#'
#' Manne, Alan. 1960. "On the Job-Shop Scheduling Problem." Operations Research 8 (2): 219-23. \doi{10.1287/opre.8.2.219}.
#'
#' Moore, Andrew, and C. G. Atkeson. 1993. "Prioritized Sweeping: Reinforcement Learning with Less Data and Less Real Time." Machine Learning 13 (1): 103–30. \doi{10.1007/BF00993104}.
#'
#' Puterman, Martin L., and Moon Chirl Shin. 1978. "Modified Policy Iteration Algorithms for Discounted Markov Decision Problems." Management Science 24: 1127-37. \doi{10.1287/mnsc.24.11.1127}.
#'
#' Rummery, G., and Mahesan Niranjan. 1994. "On-Line Q-Learning Using Connectionist Systems." Techreport CUED/F-INFENG/TR 166. Cambridge University Engineering Department.
#'
#' Russell, Stuart J., and Peter Norvig. 2020. Artificial Intelligence: A Modern Approach (4th Edition). Pearson. [http://aima.cs.berkeley.edu/](http://aima.cs.berkeley.edu/).
#'
#' Sutton, R. 1988. "Learning to Predict by the Method of Temporal Differences." Machine Learning 3: 9-44. [https://link.springer.com/article/10.1007/BF00115009](https://link.springer.com/article/10.1007/BF00115009).
#'
#' Sutton, Richard S., and Andrew G. Barto. 2018. Reinforcement Learning: An Introduction. Second. The MIT Press. [http://incompleteideas.net/book/the-book-2nd.html](http://incompleteideas.net/book/the-book-2nd.html).
#'
#' Watkins, Christopher J. C. H., and Peter Dayan. 1992. "Q-Learning." Machine Learning 8 (3): 279-92. \doi{10.1007/BF00992698}.
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
#' gw_plot(maze_solved)
#'
#' # Use linear programming
#' maze_solved <- solve_MDP(Maze, method = "lp")
#' maze_solved
#' policy(maze_solved)
#'
#' # use prioritized sweeping (which is known to be fast for mazes)
#' maze_solved <- solve_MDP(Maze, method = "prioritized_sweeping")
#' policy(maze_solved)
#'
#' # finite horizon
#' maze_solved <- solve_MDP(Maze, method = "value_iteration", horizon = 3)
#' policy(maze_solved)
#' gw_plot(maze_solved, epoch = 1)
#' gw_plot(maze_solved, epoch = 2)
#' gw_plot(maze_solved, epoch = 3)
#'
#' # create a random policy where action n is very likely and approximate
#' #  the value function. We change the discount factor to .9 for this.
#' Maze_discounted <- Maze
#' Maze_discounted$discount <- .9
#' pi <- random_policy(Maze_discounted,
#'                     prob = c(n = .7, e = .1, s = .1, w = 0.1)
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
#' maze_learned <- solve_MDP(Maze, method = "q_learning", n = 100, horizon = 100)
#' maze_learned
#'
#' maze_learned$solution
#' policy(maze_learned)
#' plot_value_function(maze_learned)
#' gw_plot(maze_learned)
#' @export
solve_MDP <- function(model, method = "value_iteration", ...) {
  if (!inherits(model, "MDP")) {
    stop("x needs to be an MDP!")
  }
  
  methods_DP <- c(
    "value_iteration",
    "policy_iteration",
    "prioritized_sweeping",
    "GenPS",
    "PS_error",
    "PS_random"
  )
  methods_LP <- c("lp")
  methods_TD <- c("sarsa", "q_learning", "expected_sarsa")
  methods_MC <- c("MC_exploring_starts", "MC_on_policy", "MC_off_policy")
  methods_sampling <- c("q_planning")
  method <- match.arg(method,
                      c(
                        methods_DP,
                        methods_LP,
                        methods_TD,
                        methods_MC,
                        methods_sampling
                      ))
  
  if (method %in% methods_DP) {
    return(solve_MDP_DP(model, method, ...))
  }
  
  if (method %in% methods_LP) {
    return(solve_MDP_LP(model, method, ...))
  }
  
  if (method %in% methods_TD) {
    return(solve_MDP_TD(model, method, ...))
  }
  
  if (method %in% methods_MC) {
    return(solve_MDP_MC(model, method, ...))
  }
  
  if (method %in% methods_sampling) {
    return(solve_MDP_sampling(model, method, ...))
  }
  
  # we should not get here!
  stop("Unknown method!")
}
