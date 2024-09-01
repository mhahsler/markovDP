## context("solve_MDP")

data("Maze")

methods_DP <- c("value_iteration", "policy_iteration", "prioritized_sweeping")
methods_LP <- c("lp")
methods_TD <- c("sarsa", "q_learning", "expected_sarsa")
methods_MC <- c("MC_exploring_starts", "MC_on_policy", "MC_off_policy")
methods_sampling <- c("q_planning")

methods <- c(
  methods_DP,
  # methods_LP, ### to slow
  methods_TD,
  methods_sampling
)

for (m in methods) {
  sol <- solve_MDP(Maze, method = m)
  pol <- policy(sol)
  expect_identical(dim(pol), c(length(Maze$states), 3L))

  # check_and_fix_MDP(sol)
}

### MC methods need a max horizon and are slow! We use a low N
for (m in methods_MC) {
  sol <- solve_MDP(Maze, method = m, N = 10, horizon = 100)
  pol <- policy(sol)
  expect_identical(dim(pol), c(length(Maze$states), 3L))

  # check_and_fix_MDP(sol)
}
