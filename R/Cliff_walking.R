#' Cliff Walking Gridworld MDP
#'
#' The cliff walking gridworld MDP example from Chapter 6 of the textbook
#' "Reinforcement Learning: An Introduction."
#'
#' The cliff walking gridworld has the following layout:
#'
#' ![Cliff Walking Gridworld](cliff-walking-gridworld.png "Cliff Walking Gridworld.")
#'
#' The gridworld is represented as a 4 x 12 matrix of states.
#' The states are labeled with their x and y coordinates.
#' The start state is in the bottom left corner.
#' Each action has a reward of -1, falling off the cliff has a reward of -100 and
#' returns the agent back to the start. The episode is finished once the agent
#' reaches the absorbing goal state in the bottom right corner.
#' No discounting is used (i.e., \eqn{\gamma = 1}).
#'
#' @docType data
#' @name Cliff_walking
#' @aliases Cliff_walking cliff_walking
#' @format An object of class [MDP].
#' @keywords datasets
#' @family MDP_examples
#' @family gridworld
#' @references
#' Richard S. Sutton and Andrew G. Barto (2018). Reinforcement Learning: An Introduction
#' Second Edition, MIT Press, Cambridge, MA.
#' @examples
#' data(Cliff_walking)
#' Cliff_walking
#'
#' gw_matrix(Cliff_walking)
#' gw_matrix(Cliff_walking, what = "labels")
#'
#' # The Goal is an absorbing state
#' absorbing_states(Cliff_walking, sparse = "states")
#'
#' # visualize the transition graph
#' gw_plot_transition_graph(Cliff_walking)
#'
#' # solve using different methods
#' sol <- solve_MDP(Cliff_walking)
#' sol
#' policy(sol)
#' gw_plot(sol)
NULL
