
for (m in models) {
  gridworld_matrix(m, what = "states")
  gridworld_matrix(m, what = "index")
  gridworld_matrix(m, what = "labels")
  expect_error(gridworld_matrix(m, what = "values"))
  expect_error(gridworld_matrix(m, what = "actions"))
  gridworld_matrix(m, what = "absorbing")
  gridworld_matrix(m, what = "unreachable")
}
