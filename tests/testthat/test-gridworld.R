models <- c(models_matrix, models_trans_function)

whats <- c("states", "index", "labels", 
           "absorbing", "unreachable")

for (w in whats) {
  res <- gw_matrix(models[[1]], what = w)
  for (m in models)
    expect_equal(gw_matrix(m, what = w), res)
}

# throw error for unsolved models
whats <- c("values", "actions") 
for (w in whats) {
  for (m in models)
    expect_error(gw_matrix(m, what = w))
}

