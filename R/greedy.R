#' Greedy Actions and Policies
#'
#' Extract a greedy policy or select a greedy action 
#' from a solved model or a Q matrix.
#'
#' @family MDP
#' @family policy
#' @author Michael Hahsler
#'
#' @param x a solved MDP model or a Q matrix.
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
#' Q <- Q_values(Maze, V)
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
#' 
#' @param epsilon an `epsilon > 0` applies an epsilon-greedy policy.
#' @param prob logical; return a probability distribution over the actions.
#' @return `greedy_action()` returns the action with the highest q-value
#'    for state `s`. If `prob = TRUE`, then a vector with
#'    the probability for each action is returned.
#' @export
greedy_action <- function(x,  
                          s,
                          epsilon = 0,
                          prob = FALSE) {
  UseMethod("greedy_action")
}

#' @export
greedy_action.matrix <-
  function(x,
           s,
           epsilon = 0,
           prob = FALSE) {
    Q <- x
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

#' @export
greedy_action.MDP <-
  function(x,
           s,
           epsilon = 0,
           prob = FALSE) {
    greedy_action(Q_values(x), s, epsilon, prob)
  }


#' @rdname greedy_action
#' @return `greedy_policy()` returns the greedy policy given `Q`.
#' @export
greedy_policy <- function(x) {  
  UseMethod("greedy_policy")
}

#' @export
greedy_policy.matrix <-
  function(x) {
    Q <- x
    A <- colnames(Q)
    a <- apply(Q, MARGIN = 1, which.max.random)
    
    data.frame(
      state = rownames(Q),
      V = Q[cbind(seq_len(nrow(Q)), a)],
      action = factor(A[a], levels = A),
      row.names = NULL
    )
  }

#' @export
greedy_policy.MDP <-
  function(x) {
    greedy_policy(Q_values(x))
  }