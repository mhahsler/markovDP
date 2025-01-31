#' Bellman Update and Bellman operator
#'
#' Update the value function with a Bellman update.
#' 
#' The Bellman update updates a value function given the model 
#' by applying the Bellman equation as an update rule for each state:
#'
#' \deqn{v_{k+1}(s) \leftarrow \max_{a \in \mathcal{A}(s)} \sum_{s'} p(s' | s,a) [r(s,a, s') + \gamma v_k(s')]}
#'
#' The Bellman update moves the estimated value function \eqn{V} closer to the 
#' optimal value function \eqn{v_*}.
#'
#' The Bellman operator \eqn{B_\pi} updates a value function given the model, 
#' and a policy \eqn{\pi}:
#'
#' \deqn{(B_\pi v)(s) = \sum_{a \in \mathcal{A}} \pi(a|s) \sum_{s'} p(s' | s,a) [r(s,a,s') + \gamma v(s')]}
#'
#' The Bellman error is \eqn{\delta = B_\pi v - v}.
#' The Bellman operator reduces the Bellman error and moves the value function 
#' closer to the fixed point of the true value function:
#'  
#' \deqn{v_\pi = B_\pi v_\pi.}
#'
#' @family MDP
#' @family policy
#' @author Michael Hahsler
#'
#' @param model an MDP problem specification.
#' @param pi a policy as a data.frame with at least columns for states and action. If `NULL`,
#'     then the policy in model is used.
#' @param V a vector representing the value function. A single 0 can be 
#'    used as a shorthand for a value function with all 0s.
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
#' bellman_update(Maze, V = 0)
#'
#' # perform simple value iteration for 10 iterations
#' V <- 0
#' for (i in seq(10))
#'   V <- bellman_update(Maze, V)$V
#'   
#' V
#' @export
bellman_update <- function(model, V) {
  if (!inherits(model, "MDP")) {
    stop("'model' needs to be of class 'MDP'.")
  }
  
  S <- S(model)
  A <- A(model)
  gamma <- model$discount %||% 1
 
  if (length(V) == 1 && all(V == 0))
    V <- numeric(length(S))
   
  Q <- sapply(A, FUN = function(a) {
    P <- transition_matrix(model, a, sparse = FALSE)
    R <- reward_matrix(model, a, sparse = FALSE)
   
    ### TODO: add a sparse version
    #QV_R(V, P, R, gamma)
    QV_cpp(V, P, R, gamma) 
  })

  rownames(Q) <- S
   
  m <- apply(Q, MARGIN = 1, which.max.random)
  
  list(V = Q[cbind(seq_along(S), m)] , 
       pi = normalize_action(m, model),
       Q = Q)
}

# Bellman update for a single state
.bellman_state_update <- function(model, s, V = 0) {
  S <- S(model)
  A <- A(model)
  gamma <- model$discount %||% 1
  if (length(V) == 1 && all(V == 0))
    V <- numeric(length(S))
  
  Qs <- sapply(A, FUN = function(a) {
    sum(transition_matrix(model, a, s, sparse = FALSE) * 
          (reward_matrix(model, a, s, sparse = FALSE) + gamma * V));  
  })
  
  m <- which.max.random(Qs)
  
  list(V = unname(Qs[m]), 
       pi = normalize_action(m, model)
  )
}

# Calculate Bellman updated state values for a given action.

# this is for a given action (the matrices are for that action) and returns a vector
# QV_cpp is the C++ implementation

# P and R are a matrix for a given action
QV_R <- function(V, P, R, gamma) {
  sapply(seq_len(nrow(P)), FUN = function(s) {
    sum(P[s, ] * (R[s, ] + gamma * V), na.rm = TRUE)
  })
}
  
# P and R are from a model (could be functions)
# This implementation needs the least memory but is slow
QV_R_model <- function(V, a, model, gamma) {
  sapply(S(model), FUN = function(s) {
    sum(transition_matrix(model, a, s) * (reward_matrix(model, a, s) + gamma * V))
  })
}


#' @rdname bellman_update
#' @export
bellman_operator <- function(model, pi, V) {
  p_pi <- induced_transition_matrix(model, pi)
  r_pi <- induced_reward_matrix(model, pi)
  r_pi <- rowSums(p_pi * r_pi)
  
  r_pi + gamma * p_pi %*% V
}

