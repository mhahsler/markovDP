# test MDP

m <- gw_maze_MDP(c(3, 3), start = "s(1,1)", goal = "s(3,3)")

# construct state features as the x/y coordinates in the gridworld
state_features <- gw_s2rc(S(m))
state_features

m <- add_linear_approx_Q_function(m, state_features)

# constructed state-action features (X) and approximate Q function
# and gradient
m$approx_Q_function

warning("Fix this!")

sol <- solve_MDP_APPROX(m, horizon = 100, n = 10,
                     alpha = 0.01, epsilon = .05)

gw_plot(sol)

sample_MDP(sol, n = 10, horizon = 100, verbose = interactive())

# test MDPE

m <- gw_maze_MDPE(c(3, 3), start = s(1,1), goal = s(3,3))
m <- add_linear_approx_Q_function(m)
m$approx_Q_function

sol <- solve_MDP_APPROX(m, horizon = 100, n = 10,
                     alpha = 0.01, epsilon = .05)

gw_plot(sol)

sample_MDP(sol, n = 10, horizon = 100, verbose = interactive())
