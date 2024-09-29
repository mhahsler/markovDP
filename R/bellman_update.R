#' Bellman Update
#'
#' Update the value function with a Bellman update.
#' 
#' The Bellman updates a value function given the model defining
#' \eqn{T}, \eqn{\gamma} and \eqn{R}  
#' by applying the Bellman equation as an update rule for each state:
#'
#' \deqn{U_{k+1}(s) \leftarrow \max_{a \in A(s)} \sum_{s'} T(s' | s,a) [R(s,a) + \gamma U_k(s')]}
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
#'    
#' @return an updates vector with state values (U).
#'
#' @references
#' Sutton, R. S., Barto, A. G. (2020). Reinforcement Learning: An Introduction.
#' Second edition. The MIT Press.
#'
#' @examples
#' data(Maze)
#' Maze
#'
#' # single Bellman update from a all 0 value function
#' bellman_update(Maze, U = 0)
#'
#' # perform simple value iteration 
#' U <- 0
#' for (i in seq(10))
#'   U <- bellman_update(Maze, U)
#'   
#' U
#' @export
bellman_update <- function(model, U) {
  S <- model$states
  A <- model$actions
  GAMMA <- model$discount %||% 1
  
  Qs <- outer(S,
              A,
              .QV_func_vec,
              model,
              GAMMA,
              U)
  
  m <- apply(Qs, MARGIN = 1, which.max.random)
  Qs[cbind(seq_along(S), m)]
}
