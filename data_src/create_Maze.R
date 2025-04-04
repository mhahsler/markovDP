library(markovDP)
# The problem can be loaded using data(Maze).

# Here is the complete problem definition:
# the wall at s(2,2) is blocked
gw <- gw_init(dim = c(3, 4),
                     start = "s(3,1)",
                     goal = "s(1,4)",
                     absorbing_states = c("s(1,4)", "s(2,4)"),
                     blocked_states = "s(2,2)",
                     state_labels = list(
                         "s(3,1)" = "Start",
                         "s(2,4)" = "-1",
                         "s(1,4)" = "Goal: +1"
                         )
                     )
gw_matrix(gw)
gw_matrix(gw, what = "labels")

# gw_init has created the following information
str(gw)

# the transition function is stochastic so we cannot use the standard
# gridworld gw$transition_prob() function and have to replace it
T <- function(model, action, start.state) {
  action <- match.arg(action, choices = model$actions)
  
  P <- structure(numeric(length(model$states)), names = model$states)
  
  # absorbing states
  if (start.state %in% model$info$absorbing_states) {
    P[start.state] <- 1
    return(P)
  }
  
  if (action %in% c("up", "down")) {
    error_direction <- c("right", "left")
  } else {
    error_direction <- c("up", "down")
  }
  
  rc <- gw_s2rc(start.state)
  delta <- list(
    up = c(-1, 0),
    down = c(+1, 0),
    right = c(0, +1),
    left = c(0, -1)
  )
  
  # there are 3 directions. For blocked directions, stay in place
  # 1) action works .8
  rc_new <- gw_rc2s(rc + delta[[action]])
  if (rc_new %in% model$states)
    P[rc_new] <- .8
  else
    P[start.state] <- .8
  
  # 2) off to the right .1
  rc_new <- gw_rc2s(rc + delta[[error_direction[1]]])
  if (rc_new %in% model$states)
    P[rc_new] <- .1
  else
    P[start.state] <-  P[start.state] + .1
  
  # 3) off to the left .1
  rc_new <- gw_rc2s(rc + delta[[error_direction[2]]])
  if (rc_new %in% model$states)
    P[rc_new] <- .1
  else
    P[start.state] <-  P[start.state] + .1
  
  P
}

T(gw, "up", "s(3,1)")


R <- rbind(
  R_(                         value = -0.04),
  R_(end.state = "s(2,4)",    value = -1-0.04),
  R_(end.state = "s(1,4)",    value = +1-0.04),
  R_(start.state = "s(2,4)",  value = 0),
  R_(start.state = "s(1,4)",  value = 0)
)


Maze <- MDP(
  name = "Stuart Russell's 3x4 Maze",
  discount = 1,
  horizon = Inf,
  states = gw$states,
  actions = gw$actions,
  start = "s(3,1)",
  transition_prob = T,
  reward = R,
  info = gw$info
)

Maze <- normalize_MDP(Maze)

Maze

save(Maze, file = "data/Maze.rda")

