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
#' Maze_with_policy <- add_policy(Maze, pi_greedy)
#' gw_plot(Maze_with_policy, main = "Maze: Greedy Policy")
#'
#' # find the greedy/ epsilon-greedy action for the top-left corner state 
#' greedy_action(Maze, "s(1,1)", Q, epsilon = 0, prob = FALSE)
#' greedy_action(Maze, "s(1,1)", Q, epsilon = 0, prob = TRUE)
#' greedy_action(Maze, "s(1,1)", Q, epsilon = .1, prob = TRUE)
#'
#' # we can also specify a model with a policy and use the internal Q-values 
#' greedy_action(Maze_with_policy, "s(1,1)", epsilon = .1, prob = TRUE)
#' 
#' @param Q an optional Q-matrix.
#' @param epsilon an `epsilon > 0` applies an epsilon-greedy policy.
#' @param prob logical; return a probability distribution over the actions.
#' @return 
#'    * `greedy_action()` returns the action with the highest q-value
#'    for state `s`. If `prob = TRUE`, then a vector with
#'    the probability for each action is returned.
#'    * `greedy_policy()` returns a data.frame with the policy.
#' @export
greedy_action <- function(x,  
                          s,
                          Q = NULL,
                          epsilon = 0,
                          prob = FALSE) {
  UseMethod("greedy_action")
}

# For a Q-matrix
greedy_action_int <-
  function(x,
           s,
           epsilon = 0,
           prob = FALSE) {
    if (missing(s))
      stop("state s needs to be specified!")
    
    Q <- x
    
    A <- colnames(Q)
    if (is.null(A))
      stop("Q matrix needs actions as column names.")
    
    # R = -Inf means unavailable action
    available_A <- which(Q[s, ] != -Inf)
    
    if (!prob) {
      if (epsilon == 0 ||
          length(available_A) == 1L || runif(1) > epsilon) {
        a <- available_A[which.max.random(Q[s, available_A])]
      } else {
        a <- sample(available_A, size = 1L)
      }
     
      a <- factor(a, levels = seq_along(A), labels = A) 
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
           Q = NULL,
           epsilon = 0,
           prob = FALSE) {
    Q <- Q %||% Q_values(x)
    greedy_action_int(Q, normalize_state_id(s, x), epsilon, prob)
  }

#' @export
greedy_action.MDPTF <- 
  function(x,
           s,
           Q = NULL,
           epsilon = 0,
           prob = FALSE) {
    Q <- Q %||% rbind(approx_Q_value(x, s))
    greedy_action_int(Q, 1L, epsilon, prob)
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

#' @export
greedy_policy.MDPTF <-
  function(x) {
    if (!is.null(S(x)))
      stop("MDPTF does not specify the state space!")
    
    approx_greedy_policy(Q_values(x))
  }