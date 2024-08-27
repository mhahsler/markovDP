#' Q-Values and Greedy Policies
#'
#' Implementation several functions useful to deal with Q-values for MDPs
#' which maps a state/action pair to a utility value.
#'
#' Implemented functions are:
#'
#' * `q_values()` calculates (approximates)
#'   Q-values for a given model and value function using the Bellman
#'   optimality equation:
#'
#'   \deqn{q(s,a) = \sum_{s'} T(s'|s,a) [R(s,a) + \gamma U(s')]}
#'
#'   Q-values are calculated if \eqn{U = U^*}, the optimal value function
#'   otherwise we get an approximation.
#'   Q-values can be used as the input for several other functions.
#'
#' * `greedy_action()` returns the action with the largest Q-value given a
#'    state.
#'
#' * `greedy_policy()`
#'    generates a greedy policy using Q-values.
#'
#' @name q_values
#' @aliases q_values
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
#' @param Q an action value function with Q-values as a state by action matrix.
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
#' u <- policy_evaluation(Maze, pi_random)
#' q <- q_values(Maze, U = u)
#'
#' # get the greedy policy form the q-values
#' pi_greedy <- greedy_policy(q)
#' pi_greedy
#' gridworld_plot(add_policy(Maze, pi_greedy), main = "Maze: Greedy Policy")
#'
#' greedy_action(q, "s(3,1)", epsilon = 0, prob = FALSE)
#' greedy_action(q, "s(3,1)", epsilon = 0, prob = TRUE)
#' greedy_action(q, "s(3,1)", epsilon = .1, prob = TRUE)
NULL

# Calculate Q-Function from U
# the (optimal) state-action value function Q_s(a,k) is the expected total reward
# from stage k onward, if we choose a_k = a and then proceed optimally (given by U).
.QV <-
  function(s, a, P, R, GAMMA, U) {
    sum(P[[a]][s, ] * (R[[a]][s, ] + GAMMA * U), 
        na.rm = TRUE)
  }
.QV_vec <- Vectorize(.QV, vectorize.args = c("s", "a"))


#' @rdname q_values
#' @return `q_values()` returns a state by action matrix specifying the Q-function,
#'   i.e., the action value for executing each action in each state. The Q-values
#'   are calculated from the value function (U) and the transition model.
#' @export
q_values <- function(model, U = NULL) {
  if (!inherits(model, "MDP")) {
    stop("'model' needs to be of class 'MDP'.")
  }

  S <- model$states
  A <- model$actions
  P <- transition_matrix(model, sparse = TRUE)
  # TODO: sparse = TRUE returns a dataframe. This uses a lot of memory!
  R <- reward_matrix(model, sparse = FALSE)
  policy <- model$solution$policy[[1]]
  GAMMA <- model$discount

  if (is.null(U)) {
    ## Q is stored in model
    if (!is.null(model$solution$Q))
      return(model$solution$Q)
    
    ## U is stored in model
    if (!is.null(policy) || any(is.na(policy$U))) {
      U <- policy$U
    } else {
      stop("'model' does not contain state utilities (it is unsolved). You need to specify U.")
    }
  }
  
  structure(outer(S, A, .QV_vec, P, R, GAMMA, U), dimnames = list(S, A))
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
    data.frame(
      state = rownames(Q),
      U = apply(Q, MARGIN = 1, max),
      action = factor(A[apply(Q, MARGIN = 1, which.max.random)], levels = A),
      row.names = NULL
    )
  }
