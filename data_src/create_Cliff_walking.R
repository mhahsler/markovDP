library(markovDP)

# Here is the complete problem definition:
NROWS <- 4
NCOLS <- 12

START <- "s(4,1)"
GOAL <- "s(4,12)"

gw <- gw_init(dim = c(NROWS, NCOLS), 
                     start = START,
                     goal = GOAL,
                     blocked_states =
                       c("s(4,2)", "s(4,3)", "s(4,4)", "s(4,5)", 
                         "s(4,6)", "s(4,7)", "s(4,8)", "s(4,9)", 
                         "s(4,10)", "s(4,11)"))
gw_matrix(gw)


T <- function(model, action, start.state) {
  P <- structure(numeric(length(model$states)), names = model$states)
  action <- match.arg(action, choices = model$actions)
  
  # GOAL is absorbing
  if (start.state == model$info$goal) {
    P[start.state] <- 1
    return(P)
  }
  
  # the cliff leads to START instead
  start_rc <- gw_s2rc(start.state)
  if ((action == "down" &&
       start_rc[1] == 3 &&
       start_rc[2] >= 2 && start_rc[2] <= 11) ||
      (start.state == "s(4,1)" && action == "right") ||
      (start.state == "s(4,12)" && action == "left")) {
    P[model$info$start] <- 1
    return(P)
  }
  
  # rest of the actions are normal
  return(gw_transition_prob(model, action, start.state))
}

T(gw, "up", "s(4,1)")
T(gw, "down", "s(3,2)")
T(gw, "down", GOAL)

R <- rbind(
  R_(value = -1),
  R_(end.state = START, value = -100),
  R_(
    action = "down",
    start.state = "s(3,1)",
    end.state = START,
    value = -1
  ),
  R_(
    start.state = GOAL,
    end.state = GOAL,
    value = 0
  )
)

R

Cliff_walking <- MDP(
  name = "Cliff Walking Gridworld",
  discount = 1,
  horizon = Inf,
  states = gw$states,
  actions = gw$actions,
  start = gw$start,
  transition_prob = T,
  reward = R,
  info = gw$info
)

Cliff_walking <- normalize_MDP(Cliff_walking)

Cliff_walking

gw_plot_transition_graph(Cliff_walking)

#Cliff_walking <- remove_unreachable_states(Cliff_walking)
#Cliff_walking

save(Cliff_walking, file = "data/Cliff_walking.rda")
