#' Solve an MDP Problem
#'
#' Implementation of value iteration, modified policy iteration and other
#' methods based on reinforcement learning techniques to solve finite
#' state space MDPs.
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
#' A list of available solvers can be found in the See Also section under 
#' "Other solvers".
#' 
#' While [MDP] model contain an explicit specification of the state space,
#' the transition probabilities and the reward structure, [MDPE] only contains a 
#' transition function. This means that only a small subset of solvers can be 
#' used for MDPEs. This currently includes only includes 
#' the solvers in [`solve_MDP_APPROX()`].  
#' 
#' @family solver
#' @family MDP
#' @family MDPE
#'
#' @param model an MDP problem specification.
#' @param method string; Composed of the algorithm family abbreviation and the algorithm
#'  separated by `:`. The algorithm families can be found in the See Also
#'  section under "Other solvers". The family abbreviation follows `solve_MDP_` in the function name. 
#' @param horizon an integer with the number of epochs for problems with a
#'   finite planning horizon. If set to `Inf`, the algorithm continues
#'   running iterations till it converges to the infinite horizon solution. If
#'   `NULL`, then the horizon specified in `model` will be used.
#' @param discount discount factor in range \eqn{(0, 1]}. If `NULL`, then the
#'   discount factor specified in `model` will be used.
#' @param continue logical; Continue with an unconverged solution specified in `model`.
#' @param matrix logical; if `TRUE` then matrices for the transition model and
#'    the reward function are taken from the model first. This can be slow if functions
#'    need to be converted or do not fit into memory if the models are large. If these
#'    components are already matrices, then this is very fast. For `FALSE`, the
#'    transition probabilities and the reward is extracted when needed. This is slower,
#'    but removes the time and memory requirements needed to calculate the matrices.
#' @param continue logical; show a progress bar with estimated time for completion.
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
#'
#' Russell, Stuart J., and Peter Norvig. 2020. Artificial Intelligence: A Modern Approach (4th Edition). Pearson. [http://aima.cs.berkeley.edu/](http://aima.cs.berkeley.edu/).
#'
#' Sutton, Richard S., and Andrew G. Barto. 2018. Reinforcement Learning: An Introduction. Second. The MIT Press. [http://incompleteideas.net/book/the-book-2nd.html](http://incompleteideas.net/book/the-book-2nd.html).
#'
#' @examples
#' data(Maze)
#' Maze
#'
#' # default is value iteration (VI)
#' maze_solved <- solve_MDP(Maze)
#' maze_solved
#' policy(maze_solved)
#'
#' # plot the value function U
#' plot_value_function(maze_solved)
#'
#' # Gridworld solutions can be visualized
#' gw_plot(maze_solved)
#'
#' @export
solve_MDP <- function(model, ...) {
  UseMethod("solve_MDP")
}

#' @rdname solve_MDP
#' @export
solve_MDP.MDP <- function(model,
                      method = "DP:VI",
                      horizon = NULL,
                      discount = NULL,
                      ...,
                      matrix = TRUE,
                      continue = FALSE,
                      verbose = FALSE,
                      progress = !verbose) {
  # remove prefix
  splt <- strsplit(method, ":")[[1]]
  prefix <- splt[1]
  method <- splt[2]
   
  func <- get(paste0("solve_MDP_", prefix))
  
  func(
    model,
    method,
    horizon,
    discount,
    ...,
    continue = continue,
    verbose = verbose,
    progress = progress
  )
}

#' @rdname solve_MDP
#' @export
solve_MDP.MDPE <- function(model,
                      method = "APPROX:semi_gradient_sarsa",
                      horizon = NULL,
                      discount = NULL,
                      ...,
                      matrix = TRUE,
                      continue = FALSE,
                      verbose = FALSE,
                      progress = !verbose) {
  # remove prefix
  splt <- strsplit(method, ":")[[1]]
  prefix <- splt[1]
  method <- splt[2]
  
  func <- get(paste0("solve_MDP_", prefix))
  
  func(
    model,
    method,
    horizon,
    discount,
    ...,
    continue = continue,
    verbose = verbose,
    progress = progress
  )
}


.prep_model <- function(model,
                        horizon = NULL,
                        discount = NULL,
                        matrix = FALSE,
                        verbose = FALSE,
                        progress = TRUE,
                        allow_MDPE = FALSE) {
  
  if (!inherits(model, "MDP") && !(allow_MDPE && inherits(model, "MDPE"))) {
    stop("'model' needs to be of class 'MDP'.")
  }
  
  model$horizon <- horizon %||% model$horizon %||% Inf
  model$discount <- discount %||% model$discount %||% 1
  
  if (inherits(model, "MDP") && matrix) {
    if (verbose)
      cat("Precomputing matrices for R and T ...")
    model <- normalize_MDP(
      model,
      sparse = NULL,
      precompute_absorbing = FALSE,
      progress = progress
    )
    
    if (verbose)
      cat(" done.\n")
  }
  
  model
}


