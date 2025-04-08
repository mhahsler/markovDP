m <- gw_maze_MDP(c(5, 5), start = "s(1,1)", goal = "s(5,5)")

benchmark <- solve_MDP(m)
gw_plot(benchmark)
policy(benchmark)
gw_matrix(benchmark, what = "value")

# goal state action is undefined!
expect_true(all(gw_matrix(benchmark, what = "action")[-25] %in% c("down", "right")))



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

set.seed(2000)
sol <- solve_MDP_APPROX(m, horizon = 100, n = 100)

policy(sol)
gw_matrix(sol, what = "value")
gw_plot(sol)

expect_true(all (gw_matrix(sol, what = "action") %in% c("down", "right")))

sample_MDP(sol, n = 10, horizon = 100, verbose = interactive())

approx_V_plot(sol)

## test MDPTF with state space
m <- gw_maze_MDPTF(c(5, 5), start = s(1,1), goal = s(5,5))

# FIXME: large n (1000) lead to anchoring problems for linear basis!
m <- add_linear_approx_Q_function(m, transformation = transformation_linear_basis)

set.seed(2000)
sol <- solve_MDP_APPROX(m, horizon = 100, n = 100)
gw_matrix(sol, what = "value")
gw_plot(sol)
policy(sol)
expect_true(all (gw_matrix(sol, what = "action") %in% c("down", "right")))

approx_V_plot(sol)

###
  
m <- add_linear_approx_Q_function(m, transformation = transformation_polynomial_basis, order = 1)
set.seed(2000)
sol <- solve_MDP_APPROX(m, horizon = 100, n = 100)
gw_matrix(sol, what = "value")
gw_plot(sol)
policy(sol)
expect_true(all (gw_matrix(sol, what = "action") %in% c("down", "right")))
  
approx_V_plot(sol)
  
###

m <- add_linear_approx_Q_function(m, transformation = transformation_RBF_basis, n = 3)
set.seed(2000)
sol <- solve_MDP_APPROX(m, horizon = 100, n = 100)

gw_matrix(sol, what = "value")
gw_plot(sol)
policy(sol)
expect_true(all (gw_matrix(sol, what = "action") %in% c("down", "right")))

approx_V_plot(sol)

##

m <- add_linear_approx_Q_function(m, transformation = transformation_fourier_basis, order = 1)
set.seed(2000)
sol <- solve_MDP_APPROX(m, horizon = 100, n = 100)
gw_matrix(sol, what = "value")
gw_plot(sol)
policy(sol)
expect_true(all (gw_matrix(sol, what = "action") %in% c("down", "right")))

approx_V_plot(sol)

# test approx without state space
states <- m$states
m$states <- NULL
S(m)

# FIXME: large n (1000) lead to anchoring problems for linear basis!
m <- add_linear_approx_Q_function(m, transformation = transformation_linear_basis)
set.seed(2000)
sol <- solve_MDP_APPROX(m, horizon = 100, n = 100)
approx_V_plot(sol, 0, 5)


m <- add_linear_approx_Q_function(m, transformation = transformation_polynomial_basis, order = 1)
set.seed(2000)
sol <- solve_MDP_APPROX(m, horizon = 100, n = 100)
approx_V_plot(sol, 0, 5)
  
  
## FIXME: Better default location of RBFs?  
m <- add_linear_approx_Q_function(m, transformation = transformation_RBF_basis, n = 4)
set.seed(2000)
sol <- solve_MDP_APPROX(m, horizon = 100, n = 100)
approx_V_plot(sol, 0, 5)

##
m <- add_linear_approx_Q_function(m, transformation = transformation_fourier_basis, order = 1)
set.seed(2000)
sol <- solve_MDP_APPROX(m, horizon = 100, n = 100)
approx_V_plot(sol, 0, 5)

# cleanup
unlink("Rplots.pdf")

