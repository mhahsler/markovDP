
set.seed(1000)

m <- gw_maze_MDP(c(5, 5), start = "s(1,1)", goal = "s(5,5)")

benchmark <- solve_MDP(m)
gw_plot(benchmark)
policy(benchmark)
gw_matrix(benchmark, what = "value")

# goal state action is undefined!
expect_true(all(gw_matrix(benchmark, what = "action")[-25] %in% c("down", "right")))

m <- add_linear_approx_Q_function(m)

sol <- solve_MDP_APPROX_LAMBDA(m, method = "sarsa", horizon = 100, n = 100)


policy(sol)
gw_matrix(sol, what = "value")
gw_plot(sol)

expect_true(all (gw_matrix(sol, what = "action") %in% c("down", "right")))
