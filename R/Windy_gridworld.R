#' Windy Gridworld MDP
#' Windy Gridworld MDP
#'
#' The Windy gridworld MDP example from Chapter 6 of the textbook
#' "Reinforcement Learning: An Introduction."
#'
#' The gridworld has the following layout:
#'
#' ![Windy Gridworld](windy-gridworld.png "Windy Gridworld.")
#'
#' The grid world is represented as a 7 x 10 matrix of states.
#' In the middle region the next states are shifted upward by wind
#' (the strength in number of squares is given below each column).
#' For example, if the agent is one cell to the right of the goal,
#' then the action left takes the agent to the cell just above the goal.
#'
#' No discounting is used (i.e., \eqn{\gamma = 1}).
#'
#' @docType data
#' @name Windy_gridworld
#' @aliases windy_gridworld
#' @format An object of class [MDP].
#' @keywords datasets
#' @family MDP_examples
#' @family gridworld
#' @references
#' Richard S. Sutton and Andrew G. Barto (2018). Reinforcement Learning: An Introduction
#' Second Edition, MIT Press, Cambridge, MA.
#' @examples
#' data(Windy_gridworld)
#' Windy_gridworld
#'
#' gw_matrix(Windy_gridworld)
#' gw_matrix(Windy_gridworld, what = "labels")
#'
#' gw_plot(Windy_gridworld)
#'
#' # The Goal is an absorbing state
#' absorbing_states(Windy_gridworld, sparse = "states")
#'
#' # visualize the transition graph
#' gw_plot_transition_graph(Windy_gridworld)
#'
#' # solve using value iteration
#' sol <- solve_MDP(Windy_gridworld)
#' sol
#' policy(sol)
#' gw_plot(sol)
NULL
