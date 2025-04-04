library(markovDP)

# Here is the complete problem definition:
nrows <- 7
ncols <- 10

START = "s(4,1)"
GOAL = "s(4,8)"

gw <- gw_init(dim = c(nrows, ncols), start = START, goal = GOAL)
gw_matrix(gw)

S <- gw$states
A <- gw$actions


T <- function(model, action, start.state, end.state) {
  action <- match.arg(action, choices = model$actions)
  
  # GOAL is absorbing
  if (start.state == model$info$goal)
    return(as.integer(end.state == start.state))
  
  # Wind towards north
  wind <- c(0, 0, 0, 1, 1, 1, 2, 2, 1, 0)
  
  rc <- gw_s2rc(start.state)
  w <- wind[rc[2]]
  
  try_move <- function(a, rc) {
    rc_new <- switch(
      a,
      "up" =     c(rc[1] - 1, rc[2]),
      "down" =   c(rc[1] + 1, rc[2]),
      "left" =   c(rc[1],     rc[2] - 1),
      "right" =  c(rc[1],     rc[2] + 1)
    )
    if (rc_new[1] >= 1 &&
        rc_new[1] <= 7 && rc_new[2] >= 1 && rc_new[2] <= 10)
      return(rc_new)
    else
      return(rc)
  }
  
  rc <- try_move(action, rc) 
    
  for (i in seq_len(w)) 
    rc <- try_move("up", rc) 
      
  new.state <- gw_rc2s(rc)
  
  
  return(as.integer(new.state == end.state))
}

T(gw, "left", "s(4,9)", "s(3,8)")
T(gw, "left", "s(4,9)", "s(3,9)")
T(gw, "right", "s(1,1)", "s(1,2)")
T(gw, "right", GOAL, GOAL)
T(gw, "up", "s(6,8)", "s(3,8)")

R <- rbind(
  R_(value = -1),
  R_(
    start.state = GOAL,
    end.state = GOAL,
    value = 0
  )
)

R

Windy_gridworld <- MDP(
  name = "Windy Gridworld",
  discount = 1,
  horizon = Inf,
  states = S,
  actions = A,
  start = START,
  transition_prob = T,
  reward = R,
  info = gw$info
)


Windy_gridworld <- normalize_MDP(Windy_gridworld)


gw_matrix(Windy_gridworld, what = "unreachable")
gw_matrix(Windy_gridworld)
gw_matrix(Windy_gridworld, what = "labels")
gw_plot_transition_graph(Windy_gridworld)
gw_plot(Windy_gridworld)

sol <- solve_MDP(Windy_gridworld)

gw_plot(sol)




save(Windy_gridworld, file = "data/Windy_gridworld.rda")
