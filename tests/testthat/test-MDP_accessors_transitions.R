# dense get all
correct <- transition_matrix(Maze_orig, sparse = FALSE)

for (m in models) {
  tr <- transition_matrix(m, sparse = FALSE)
  expect_equal(tr, correct)
}

# sparse get action matrix
correct <- transition_matrix(Maze_orig, 1, sparse = TRUE)

for (m in models) {
  tr <- transition_matrix(m, 1, sparse = TRUE)
  expect_equal(tr, correct)
}


### action/row
correct <- transition_matrix(Maze_orig, "up", 1, sparse = TRUE)

for (m in models) {
  tr <- transition_matrix(m, "up", 1, sparse = TRUE)
  expect_equal(tr, correct)
}

correct <- transition_matrix(Maze_orig, "up", 1:2, sparse = TRUE)

for (m in models) {
  tr <- transition_matrix(m, "up", 1:2, sparse = TRUE)
  expect_equal(tr, correct)
}


### action/col
correct <- transition_matrix(Maze_orig, "up", , 2, sparse = TRUE)

for (m in models) {
  tr <- transition_matrix(m, "up", , 2, sparse = TRUE)
  expect_equal(tr, correct)
}

### action/row/col
correct <- transition_matrix(Maze_orig, "up", 1, 1)
for (m in models) {
  tr <- transition_matrix(m, "up", 1, 1)
  expect_equal(tr, correct)
}


# simplify (when action is missing)
for (m in models) {
  (tr <- transition_matrix(m, NULL, 1, 1))
  expect_true(inherits(tr, "list"))
  (tr <- transition_matrix(m, NULL, 1:2, 1))
  expect_true(inherits(tr, "list"))
  (tr <- transition_matrix(m, NULL, 1, 1:2))
  expect_true(inherits(tr, "list"))
  (tr <- transition_matrix(m, NULL, 1, 1:2))
  expect_true(inherits(tr, "list"))
  (tr <- transition_matrix(m, NULL, 1:2, 1:2))
  expect_true(inherits(tr, "list"))
}

res <- transition_matrix(models[[1]], NULL, 1, 1, simplify = TRUE)
expect_true(inherits(res, "numeric"))
for (m in models)
  expect_identical(transition_matrix(m, NULL, 1, 1, simplify = TRUE), res)
  
res <- transition_matrix(models[[1]], NULL, 1:2, 1, simplify = TRUE)
expect_true(inherits(res, "matrix"))
for (m in models)
  expect_identical(transition_matrix(m, NULL, 1:2, 1, simplify = TRUE, sparse = FALSE), res)
  
res <- transition_matrix(models[[1]], NULL, 1, 1:2, simplify = TRUE)
expect_true(inherits(res, "matrix"))
for (m in models)
  expect_identical(transition_matrix(m, NULL, 1, 1:2, simplify = TRUE, sparse = FALSE), res)
  
res <- transition_matrix(models[[1]], NULL, 1:2, 1:2, simplify = TRUE)
expect_true(inherits(res, "list"))
for (m in models)
  expect_identical(transition_matrix(m, NULL, 1:2, 1:2, simplify = TRUE, 
                                     sparse = FALSE), res)


