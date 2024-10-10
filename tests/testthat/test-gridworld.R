
for (m in models) {
  gw_matrix(m, what = "states")
  gw_matrix(m, what = "index")
  gw_matrix(m, what = "labels")
  expect_error(gw_matrix(m, what = "values"))
  expect_error(gw_matrix(m, what = "actions"))
  gw_matrix(m, what = "absorbing")
  gw_matrix(m, what = "unreachable")
}
