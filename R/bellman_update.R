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
#'    (expected sum of discounted rewards from that point on). A single 0 can be 
#'    used as a shorthand for a value function with all 0s.
#' @param return_Q logical; return also the Q matrix.
#'    
#' @return a list with the updated state value vector U and the taken actions
#'    pi.
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
#' # perform simple value iteration for 10 iterations
#' U <- list(U = 0)
#' for (i in seq(10))
#'   U <- bellman_update(Maze, U$U)
#'   
#' U
#' @export
bellman_update <- function(model, U, return_Q = FALSE) {
  if (!inherits(model, "MDP")) {
    stop("'model' needs to be of class 'MDP'.")
  }
  
  S <- model$states
  A <- model$actions
  GAMMA <- model$discount %||% 1
  
 if (length(U) == 1 && all(U == 0))
   U <- numeric(length(S))
   
  Q <- sapply(A, FUN = function(a) {
    P <- transition_matrix(model, a, sparse = FALSE)
    R <- reward_matrix(model, a, sparse = FALSE)
   
    ### TODO: add a sparse version
    QV_cpp(U, P, R, GAMMA) 
  })

  rownames(Q) <- S
   
  m <- apply(Q, MARGIN = 1, which.max.random)
  
  if (!return_Q)
    list(U = Q[cbind(seq_along(S), m)] , 
         pi = factor(m, levels = seq_along(A), labels = A))
  else
    list(U = Q[cbind(seq_along(S), m)] , 
         pi = factor(m, levels = seq_along(A), labels = A),
         Q = Q)
}

# Bellman update for a single state
.bellman_state_update <- function(model, s, U = 0) {
  S <- model$states
  A <- model$actions
  GAMMA <- model$discount %||% 1
  if (length(U) == 1 && all(U == 0))
    U <- numeric(length(S))
  
  Qs <- sapply(A, FUN = function(a) {
    sum(transition_matrix(model, a, s, sparse = FALSE) * 
          (reward_matrix(model, a, s, sparse = FALSE) + GAMMA * U));  
  })
  
  m <- which.max.random(Qs)
  
  list(U = unname(Qs[m]), 
       pi = factor(m, levels = seq_along(A), labels = A))
}

# Calculate Bellman updated state values for a given action.

# this is for a given action (the matrices are for that action) and returns a vector
# QV_cpp is the C++ implementation

# P and R are a matrix for a given action
QV_R <- function(U, P, R, GAMMA) {
  sapply(seq_len(nrow(P)), FUN = function(s) {
    sum(P[s, ] * (R[s, ] + GAMMA * U))
  })
}
  
# P and R are from a model (could be functions)
# This implementation needs the least memory but is slow
QV_R_model <- function(U, a, model, GAMMA) {
  sapply(model$states, FUN = function(s) {
    sum(transition_matrix(model, a, s) * (reward_matrix(model, a, s) + GAMMA * U))
  })
}

