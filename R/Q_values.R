#' Q-Values
#'
#' Several useful functions to deal with Q-values (action values)
#' which map each state/action pair to a utility value.
#'
#' Implemented functions are:
#'
#' * `Q_values()` approximates
#'   Q values for a given model and value function using the Bellman
#'   optimality equation:
#'
#'   \deqn{q_*(s,a) = \sum_{s'} p(s'|s,a) [r(s,a,s') + \gamma v_*(s')]}
#'
#'   Exact Q values are calculated if \eqn{v = v_*}, the optimal value function,
#'   otherwise we get an approximation that might not 
#'   be consistent with \eqn{v} or the implied policy.
#'   Q values can be used as the input for several other functions.
#'
#' * `Q_zero()` and `Q_random()` create initial Q value matrices for algorithms.
#'
#' @name Q_values
#' @aliases Q_values
#'
#' @family MDP
#' @family policy
#' @author Michael Hahsler
#'
#' @param model an MDP problem specification.
#' @param V the state values. If `model` is a solved model, then the state
#'    values are taken from the solution.
#'
#' @references
#' Sutton, R. S., Barto, A. G. (2020). Reinforcement Learning: An Introduction.
#' Second edition. The MIT Press.
#'
#' @examples
#' data(Maze)
#' Maze
#'
#' # create a random policy and calculate q-values
#' pi_random <- random_policy(Maze)
#' pi_random
#' 
#' V <- policy_evaluation(Maze, pi_random)
#' V
#' 
#' # calculate Q values
#' Q <- Q_values(Maze, V)
#' Q
NULL

#' @rdname Q_values
#' @return `Q_values()` returns a state by action matrix specifying the Q-function,
#'   i.e., the action value for executing each action in each state. The Q-values
#'   are calculated from the value function (U) and the transition model.
#' @export
Q_values <- function(model, V = NULL) {
  if (is.null(V)) {
    V <- policy(model)$V
  }
  
  Q <- bellman_update(model, V)$Q
  
  Q
}

#' @rdname Q_values
#' @param value value to initialize the Q-value matrix. Default is 0.
#' @return `Q_zero()` and `Q_random` return a matrix with q-values.
#' @export
Q_zero <- function(model, value = 0) {
  S <- S(model)
  A <- A(model)
  
  matrix(value,
         nrow = length(S),
         ncol = length(A),
         dimnames = list(S, A)
  )
}

#' @rdname Q_values
#' @param min,max range of the random values
#' @export
Q_random <- function(model, min = 1e-6, max = 1) {
  S <- S(model)
  A <- A(model)
  
  matrix(runif(length(S) * length(A), min, max),
         nrow = length(S),
         ncol = length(A),
         dimnames = list(S, A)
  )
}

# internal function to initialize the Q-value matrix
# param Q:  
#   * complete matrix: matrix
#   * 0
#   * single value
#   * a range for rand
init_Q <- function(model, Q = NULL) {
  if (is.matrix(Q))
    return(Q)
  
  if (length(Q) == 1L)
    return(Q_zero(model, value = Q))
  
  if (length(Q) == 2L)
    return(Q_random(model, min = Q[1L], max = Q[2L]))
  
  else 
    stop("Illegal definition for initializeing the Q-value matrix.")
}