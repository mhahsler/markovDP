#' Q-Values and Greedy Policies
#'
#' Several useful functions to deal with Q-values (action values)
#' which map each state/action pair to a utility value.
#'
#' Implemented functions are:
#'
#' * `q_values()` approximates
#'   Q values for a given model and value function using the Bellman
#'   optimality equation:
#'
#'   \deqn{q_*(s,a) = \sum_{s'} p(s'|s,a) [r(s,a,s') + \gamma v_*(s')]}
#'
#'   Exact Q values are calculated if \eqn{V = v_*}, the optimal value function,
#'   otherwise we get an approximation that might not 
#'   be consistent with \eqn{V} or the implied policy.
#'   Q values can be used as the input for several other functions.
#'
#' * `greedy_action()` returns the action with the largest Q value given a
#'    state.
#'
#' * `greedy_policy()`
#'    generates a greedy policy using Q values.
#'
#' @name q_values
#' @aliases q_values
#'
#' @family MDP
#' @family policy
#' @author Michael Hahsler
#'
#' @param model an MDP problem specification.
#' @param V the state values.
#'    If `model` is a solved model, then the state
#'    values are taken from the solution.
#' @param Q an action-value function with \eqn{Q(s,a)} values as a 
#'    state by action matrix.
#' @param s a state.
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
#' Q <- q_values(Maze, V)
#' Q
#'
#' # get the greedy policy form the Q values
#' pi_greedy <- greedy_policy(Q)
#' pi_greedy
#' gw_plot(add_policy(Maze, pi_greedy), main = "Maze: Greedy Policy")
#'
#' # find the greedy/ epsilon-greedy action for the top-left corner state 
#' greedy_action(Q, "s(1,1)", epsilon = 0, prob = FALSE)
#' greedy_action(Q, "s(1,1)", epsilon = 0, prob = TRUE)
#' greedy_action(Q, "s(1,1)", epsilon = .1, prob = TRUE)
NULL

#' @rdname q_values
#' @return `q_values()` returns a state by action matrix specifying the Q-function,
#'   i.e., the action value for executing each action in each state. The Q-values
#'   are calculated from the value function (U) and the transition model.
#' @export
q_values <- function(model, V = NULL) {
  if (is.null(V)) {
    V <- policy(model)$V
  }
  
  Q <- bellman_update(model, V)$Q
  
  Q
}

#' @rdname q_values
#' @param epsilon an `epsilon > 0` applies an epsilon-greedy policy.
#' @param prob logical; return a probability distribution over the actions.
#' @return `greedy_action()` returns the action with the highest q-value
#'    for state `s`. If `prob = TRUE`, then a vector with
#'    the probability for each action is returned.
#' @export
greedy_action <-
  function(Q,
           s,
           epsilon = 0,
           prob = FALSE) {
    # R = -Inf means unavailable action
    available_A <- colnames(Q)[Q[s, ] != -Inf]
    
    if (!prob) {
      if (epsilon == 0 ||
          length(available_A) == 1L || runif(1) > epsilon) {
        a <- available_A[which.max.random(Q[s, available_A])]
      } else {
        a <- sample(available_A, size = 1L)
      }
     
      a <- factor(a, levels = colnames(Q)) 
      return(a)
    }
    
    # return probabilities
    p <- structure(rep(0, ncol(Q)), names = colnames(Q))
    a <- available_A[which.max.random(Q[s, available_A])]
    p[a] <- 1 - epsilon
    p[available_A] <- p[available_A] + epsilon / length(available_A)
    
    return(p)
  }

#' @rdname q_values
#' @return `greedy_policy()` returns the greedy policy given `Q`.
#' @export
greedy_policy <-
  function(Q) {
    A <- colnames(Q)
    a <- apply(Q, MARGIN = 1, which.max.random)
    
    data.frame(
      state = rownames(Q),
      V = Q[cbind(seq_len(nrow(Q)), a)],
      action = factor(A[a], levels = A),
      row.names = NULL
    )
  }
