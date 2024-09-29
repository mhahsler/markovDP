#' Steward Russell's 4x3 Maze Gridworld MDP
#'
#' The 4x3 maze is described in Chapter 17 of the textbook
#' "Artificial Intelligence: A Modern Approach" (AIMA).
#'
#' The simple maze has the following layout:
#'
#' \preformatted{
#'     1234           Transition model:
#'    ######             .8 (action direction)
#'   1#   +#              ^
#'   2# # -#              |
#'   3#S   #         .1 <-|-> .1
#'    ######
#' }
#'
#' We represent the maze states as a gridworld matrix with 3 rows and
#' 4 columns. The states are labeled `s(row, col)` representing the position in
#' the matrix.
#' The # (state `s(2,2)`) in the middle of the maze is an obstruction and not reachable.
#' Rewards are associated with transitions. The default reward (penalty) is -0.04.
#' The start state marked with `S` is `s(3,1)`.
#' Transitioning to `+` (state `s(1,4)`) gives a reward of +1.0,
#' transitioning to `-` (state `s_(2,4)`)
#' has a reward of -1.0. Both these states are absorbing
#' (i.e., terminal) states.
#'
#' Actions are movements (`up`, `right`, `down`, `left`). The actions are
#' unreliable with a .8 chance
#' to move in the correct direction and a 0.1 chance to instead to move in a
#' perpendicular direction leading to a stochastic transition model.
#'
#' Note that the problem has reachable terminal states which leads to a proper policy
#' (that is guaranteed to reach a terminal state). This means that the solution also
#' converges without discounting (`discount = 1`).
#' @name Maze
#' @aliases Maze maze
#' @family MDP_examples
#' @family gridworld
#' @docType data
#' @format An object of class [MDP].
#' @keywords datasets
#' @references Russell,9 S. J. and Norvig, P. (2020). Artificial Intelligence:
#'  A modern approach. 4rd ed.
#' @examples
#' # The problem can be loaded using data(Maze).
#'
#' # Here is the complete problem definition.
#' 
#' # We first look at the state layout
#' gridworld_matrix(gridworld_init(dim = c(3, 4)))
#' 
#' # the wall at s(2,2) is unreachable
#' gw <- gridworld_init(dim = c(3, 4),
#'         start = "s(3,1)",
#'         goal = "s(1,4)",
#'         absorbing_states = c("s(1,4)", "s(2,4)"),
#'         unreachable_states = "s(2,2)",
#'         state_labels = list(
#'             "s(3,1)" = "Start",
#'             "s(2,4)" = "-1",
#'             "s(1,4)" = "Goal: +1")
#' )
#' gridworld_matrix(gw)
#' gridworld_matrix(gw, what = "index")
#' gridworld_matrix(gw, what = "labels")
#'
#' # gridworld_init has created the following information
#' str(gw)
#'
#' # the transition function is stochastic so we cannot use the standard
#' # gridworld gw$transition_prob() function and have to replace it
#' T <- function(model, action, start.state) {
#'   action <- match.arg(action, choices = model$actions)
#'   
#'   P <- structure(numeric(length(model$states)), names = model$states)
#'   
#'   # absorbing states
#'   if (start.state %in% model$info$absorbing_states) {
#'     P[start.state] <- 1
#'     return(P)
#'   }
#'   
#'   if (action %in% c("up", "down")) {
#'     error_direction <- c("right", "left")
#'   } else {
#'     error_direction <- c("up", "down")
#'   }
#'   
#'   rc <- gridworld_s2rc(start.state)
#'   delta <- list(
#'     up = c(-1, 0),
#'     down = c(+1, 0),
#'     right = c(0, +1),
#'     left = c(0, -1)
#'   )
#'   
#'   # there are 3 directions. For blocked directions, stay in place
#'   # 1) action works .8
#'   rc_new <- gridworld_rc2s(rc + delta[[action]])
#'   if (rc_new %in% model$states)
#'     P[rc_new] <- .8
#'   else
#'     P[start.state] <- .8
#'   
#'   # 2) off to the right .1
#'   rc_new <- gridworld_rc2s(rc + delta[[error_direction[1]]])
#'   if (rc_new %in% model$states)
#'     P[rc_new] <- .1
#'   else
#'     P[start.state] <-  P[start.state] + .1
#'   
#'   # 3) off to the left .1
#'   rc_new <- gridworld_rc2s(rc + delta[[error_direction[2]]])
#'   if (rc_new %in% model$states)
#'     P[rc_new] <- .1
#'   else
#'     P[start.state] <-  P[start.state] + .1
#'   
#'   P
#'   } 
#'
#' T(gw, "up", "s(3,1)")
#'
#' R <- rbind(
#'   R_(                         value = -0.04),
#'   R_(end.state = "s(2,4)",    value = -1),
#'   R_(end.state = "s(1,4)",    value = +1),
#'   R_(start.state = "s(2,4)",  value = 0),
#'   R_(start.state = "s(1,4)",  value = 0)
#' )
#'
#'
#' Maze <- MDP(
#'   name = "Stuart Russell's 3x4 Maze",
#'   discount = 1,
#'   horizon = Inf,
#'   states = gw$states,
#'   actions = gw$actions,
#'   start = "s(3,1)",
#'   transition_prob = T,
#'   reward = R,
#'   info = gw$info
#' )
#'
#' Maze
#'
#' str(Maze)
#'
#' gridworld_matrix(Maze)
#' gridworld_matrix(Maze, what = "labels")
#' gridworld_plot(Maze)
#'
#' # find absorbing (terminal) states
#' absorbing_states(Maze)
#'
#' maze_solved <- solve_MDP(Maze)
#' policy(maze_solved)
#'
#' gridworld_matrix(maze_solved, what = "values")
#' gridworld_matrix(maze_solved, what = "actions")
#'
#' gridworld_plot(maze_solved)
NULL

## save(Maze, file = "data/Maze.rda")
