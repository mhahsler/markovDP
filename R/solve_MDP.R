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
#' Most solvers can be interrupted using Esc/CTRL-C and will return the 
#' current solution. Solving can be continued by calling `solve_MDP` with the 
#' partial solution as the model and the parameter `continue = TRUE`. This method 
#' can also be used to reduce parameters like `alpha` or `epsilon` (see Q-learning
#' in the Examples section).  
#' 
#' Next, we describe the different types of available solvers.
#'
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
#'     `'n_step_SARSA_on_policy'`, 
#'     `'n_step_SARSA_off_policy'`,
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
#' @param verbose logical or a numeric verbose level; if set to `TRUE` or `1`, the 
#'   function displays the used algorithm parameters and progress information. 
#'   Levels `>1` provide more detailed solver output in the R console.
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
#' Puterman, Martin L. 1996. Markov decision processes: discrete stochastic dynamic programming. John Wiley & Sons.
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
#' maze_learned <- solve_MDP(Maze, method = "q_learning", 
#'     epsilon = 0.2, n = 500, horizon = 100, verbose = TRUE)
#' maze_learned
#'
#' policy(maze_learned)
#' plot_value_function(maze_learned)
#' gw_plot(maze_learned)
#' 
#' # Keep on learning, but with a reduced epsilon
#' maze_learned <- solve_MDP(maze_learned, method = "q_learning",
#'     epsilon = 0.01, n = 500, horizon = 100, continue = TRUE, verbose = TRUE)
#' 
#' policy(maze_learned)
#' plot_value_function(maze_learned)
#' gw_plot(maze_learned)
#' 
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
  methods_TD_n_step <- c("n_step_sarsa_on_policy", "n_step_sarsa_off_policy")
  methods_MC <- c("MC_exploring_starts", "MC_on_policy", "MC_off_policy")
  methods_sampling <- c("q_planning")
  method <- match.arg(method,
                      c(
                        methods_DP,
                        methods_LP,
                        methods_TD,
                        methods_TD_n_step,
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
  
  if (method %in% methods_TD_n_step) {
    return(solve_MDP_TD_n_step(model, method, ...))
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
