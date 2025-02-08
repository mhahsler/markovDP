m <- gw_maze_MDP(c(3, 3), start = "s(1,1)", goal = "s(3,3)")

benchmark <- solve_MDP(m)
gw_plot(benchmark)
policy(benchmark)
gw_matrix(benchmark, what = "value")

# test MDP

# Note: Anchoring has an issue which will lead to too large values for
# s(2,3) and s(3,2) (> 100) when n is large! This happens because going out
# of bounds and then to the goal loses 1 x movement cost but then
# thinks that it will get the inflated the approximated value for the goal and 
# not the anchored value! This will eventually going out of bounds make the best
# action!
#
# Fourier basis transformation does better with this!


# construct state features as the x/y coordinates in the gridworld
#state_features <- gw_s2rc(S(m))
#state_features
#m <- add_linear_approx_Q_function(m, state_features)
m <- add_linear_approx_Q_function(m)

# constructed state-action features (X) and approximate Q function
# and gradient
m$approx_Q_function

sol <- solve_MDP_APPROX(m, horizon = 100, n = 10,
                     alpha = 0.01, epsilon = .1)


policy(sol)
gw_matrix(sol, what = "value")
gw_plot(sol)

sample_MDP(sol, n = 10, horizon = 100, verbose = interactive())

## test MDPTF
m <- gw_maze_MDPTF(c(3, 3), start = s(1,1), goal = s(3,3))

# linear approx. has issues with anchor state-action values (the goal)
# when n increases!
m <- add_linear_approx_Q_function(m)
m$approx_Q_function

sol <- solve_MDP_APPROX(m, horizon = 100, n = 10,
                     alpha = 0.01, epsilon = .1)

gw_matrix(sol, what = "value")
gw_plot(sol)
policy(sol)

sample_MDP(sol, n = 10, horizon = 100, verbose = interactive())

# test approx without states
states <- m$states
m$states <- NULL
S(m)

m <- add_linear_approx_Q_function(m)
sol <- solve_MDP_APPROX(m, horizon = 100, n = 10,
                     alpha = 0.01, epsilon = .1)

expect_error(policy(sol))
approx_greedy_action(sol, s(1,1))

expect_error(approx_greedy_policy(sol))
action(sol, s(1,1))

A <- matrix(sapply(states, FUN = function(s) approx_greedy_action(sol, s)), ncol = 3)
A

V <- matrix(apply(sapply(states, FUN = function(s) approx_Q_value(sol, s)), MARGIN = 2, max), ncol = 3)
V

# no policy!
#gw_plot(sol)
#gw_matrix(sol, what = "value")


sample_MDP(sol, n = 10, horizon = 100, verbose = interactive())

