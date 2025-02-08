#' Define an MDP as an Agent Environment
#'
#' Defines a discrete-time agent environment. The environment's MDP is defined
#' via a transition function and states and the transition probabilities are not
#' directly specified.
#'
#' Defines a discrete-time agent environment. The environment is defined
#' via a transition
#' function using state features (i.e., using a factored state representation).
#' In comparison to the [`MDP`], where transition probabilities are directly
#' specified and potentially used by solvers, the `MDPTF` has no specification
#' of transition probabilities and the potentially infinite set of states is
#' unknown.
#'
#' To represent states, a factored state representation as a **row** vector for a single
#' state or a matrix with row vectors for a set of states are used. State labels
#' are constructed in the form `s(feature1, feature2, ...)`. Conversion between the factored
#' representation and state labels is available in [`state2features()`]
#' and [`features2state()`]. Since the state set is not directly represented,
#' **state ids are cannot be used!**
#'
#' Reinforcement learning algorithms with approximation can be used to solve
#' these problems. See: [`solve_MDP_APPROX()`].
#'
#' @family MDPTF
#'
#' @param actions a character vector specifying the names of the available
#'  actions.
#' @param transition_func A transition function receiving the current state features and
#'  returning the reward and the next state features.
#' @param discount numeric; discount rate between 0 and 1.
#' @param horizon numeric; Number of epochs.
#' @param start Specifies in which state the MDP starts.
#' @param states Optional state labels.
#' @param absorbing_states a single state or a list with absorbing states.
#' @param info A list with additional information.
#' @param name a string to identify the MDP problem.
#'
#' @return The function returns an object of class MDPTF which is list with
#'   the model specification.
#' @author Michael Hahsler
#' @examples
#' # Define a simple 5x5 maze without walls
#'
#' transition_func <- function(model, state, action) {
#'   if (all(state == s(5, 5)))
#'     return(list(reward = 0, state_prime = state))
#'
#'   action <- normalize_action_label(action, model)
#'
#'   sp <- state + switch(action,
#'     up =   c( -1, 0),
#'     down = c( +1, 0),
#'     left = c(  0,-1),
#'     right = c( 0,+1)
#'   )
#'   r <- -1
#'
#'   # check bounds
#'   if (any(sp < 1) || any(sp > 5)) {
#'       sp <- state
#'    }
#'
#'   # goal
#'   if (all(sp == s(5, 5)))
#'     r <- 100
#'
#'   return(list(reward = r, state_prime = sp))
#' }
#'
#' m <- MDPTF(actions = c("up", "right", "down", "left"),
#'           transition_func,
#'           start = s(1,1),
#'           absorbing_states = s(5, 5),
#'           name = "5x5 Maze")
#' m
#'
#' act(m, s(1,1), "down")
#'
#' # Reach the goal: reward = 100
#' act(m, s(5,4), "right")
#'
#' # Illegal action: no movement and reward = -1
#' act(m, s(1,1), "up")
#'
#' # Absorbing state: no movement and 0 reward
#' act(m, s(5,5), "up")
#'
#' # Example: Solve using Linear Feature Approximation
#'
#' # the same maze can be created with this helper and then gridworld helper
#' #    functions can be used
#' m <- gw_maze_MDPTF(c(5,5), start = "s(1,1)", goal = "s(5,5)")
#'
#' m_approx <- add_linear_approx_Q_function(m)
#' sol <- solve_MDP_APPROX(m_approx, horizon = 1000, n = 10,
#'           alpha = 0.01, epsilon = 0.8, verbose = FALSE)
#' gw_plot(sol)
#'
#' approx_greedy_action(sol, s(1,1))
#' approx_greedy_action(sol, s(5,4))
#' # to use the greedy_action interface, the state label has to be used
#' greedy_action(sol, "s(1,1)")
#'
#' approx_Q_value(sol, s(1,1))
#' approx_Q_value(sol, s(5,4))
#'
#' # Example: Solve using order-2 Fourier Basis Features
#' fourier_basis_trans <- function(x) {
#'   x <- x / 10
#'   cs <- expand.grid(0:2, 0:2)
#'   apply(cs, MARGIN = 1, FUN = function(c) cos(pi * x %*% c))
#' }
#'
#' m_approx <- add_linear_approx_Q_function(m, transformation = fourier_basis_trans)
#' sol <- solve_MDP_APPROX(m_approx, horizon = 1000, n = 10,
#'           alpha = 0.01, epsilon = 0.7, verbose = FALSE)
#' gw_plot(sol)
#' sol <- solve_MDP_APPROX(sol, horizon = 1000, n = 100,
#'           alpha = 0.01, epsilon = 0.05, verbose = FALSE, continue = TRUE)
#' gw_plot(sol)
#' @export
MDPTF <- function(actions,
                  transition_func,
                  start,
                  states = NULL,
                  absorbing_states = NULL,
                  discount = .9,
                  horizon = Inf,
                  info = NULL,
                  name = NA) {
  structure(
    list(
      name = name,
      discount = discount,
      horizon = horizon,
      actions = actions,
      states = states,
      start = normalize_state_features(start),
      absorbing_states = normalize_state_features(absorbing_states),
      transition_func = transition_func,
      info = info
    ),
    class =  c("MDPTF", "MDPE")
  )
}

#' @export
print.MDPTF <- function(x, ...) {
  writeLines(paste(paste(class(x), collapse = ", "), "-", x$name))
  
  if (!is.null(x$discount)) {
    writeLines(sprintf("  Discount factor: %s", x$discount))
  }
  
  if (!is.null(x$horizon)) {
    writeLines(sprintf("  Horizon: %s epochs", x$horizon))
  }
  
  writeLines(sprintf("  Size: %d actions", length(x$actions)))
  
  writeLines(paste0("  Start: ", shorten(
    paste(features2state(x$start), collapse = ", "), n = -10L
  )))
  
  writeLines("")
  
  writeLines(strwrap(
    paste("List components:", paste(sQuote(names(
      x
    )), collapse = ", "), "\n"),
    indent = 2,
    exdent = 4
  ))
}
