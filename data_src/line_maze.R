library(markovDP)


# single start / no restart
lm <- gw_maze_MDP(dim = c(1,9), start = "s(1,2)", 
                         goal = "s(1,5)", restart = FALSE, name = "Line Maze")
lm$actions <- c("left", "right")

#plot_transition_graph(lm, layout = cbind(seq(-1, 1, length.out = length(lm$states)), 0), 
#                      rescale = FALSE, edge.arrow.size = .2)

gw_plot_transition_graph(lm, vertex.size = 70)

s <- solve_MDP(lm)
policy(s)
gw_plot(s)


# single start / restart
lm <- gw_maze_MDP(dim = c(1,9), start = "s(1,2)", 
                         goal = "s(1,5)", restart = TRUE, name = "Line Maze")
lm$actions <- c("left", "right", "restart")

plot_transition_graph(lm, edge.arrow.size = .2)

s <- solve_MDP(lm)
policy(s)
gw_plot_policy(s)

lm <- gw_maze_MDP(dim = c(1,9), start = "uniform", 
                         goal = "s(1,5)", restart = FALSE, name = "Line Maze")
lm$actions <- c("left", "right")


plot_transition_graph(lm)

plot_transition_graph(lm, edge.arrow.size = .2, layout = cbind(seq(-1, 1, length.out = length(lm$states)), 0), 
                      rescale = FALSE)

plot_transition_graph(lm, edge.arrow.size = .2)

gw_plot_transition_graph(lm, vertex.size = 80)

s <- solve_MDP(lm)
policy(s)
gw_plot_policy(s)

# uniform / restart
lm <- gw_maze_MDP(dim = c(1,9), start = "uniform", 
                         goal = "s(1,5)", restart = TRUE, name = "Line Maze")
# remove unused actions!
lm$actions <- c("left", "right", "restart")


plot_transition_graph(lm, layout = cbind(seq(-1, 1, length.out = length(lm$states)), 0), 
                      rescale = FALSE, edge.arrow.size = .2)

plot_transition_graph(lm, edge.arrow.size = .2, layout = igraph::layout.circle)


s <- solve_MDP(lm)
policy(s)
gw_plot_policy(s)


# q-learning
s <- solve_MDP(lm, method = "q", horizon = 100, N = 100)
policy(s)

s <- solve_MDP(lm, method = "sarsa", horizon = 100, N = 100)
policy(s)
s$solution

s <- solve_MDP(lm, method = "exp", horizon = 100, N = 100)
policy(s)



#Make this into a POMDP


# line_maze <- POMDP(states = lm$states, 
#                    actions = c("left", "right"),
#                    start = lm$start,
#                    transition_prob = lm$transition_prob, 
#                    reward = reward,
#                    observations = observations, 
#                    observation_prob = observation_func, 
#                    info = lm$info,
#                    name = "Line Maze"
# )


lm_POMDP <- make_partially_observable(lm, observations = observations, observation_prob = observation_func)

gw_plot_transition_graph(lm_POMDP, vertex.size = 100)

sol <- solve_POMDP(lm_POMDP)
