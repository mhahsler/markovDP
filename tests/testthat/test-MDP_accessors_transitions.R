# dense get all
(tr <- transition_matrix(Maze_orig, sparse = FALSE))
(tr_dense <- transition_matrix(Maze_dense, sparse = FALSE))
(tr_sparse <- transition_matrix(Maze_sparse, sparse = FALSE))
(tr_function2 <- transition_matrix(Maze_function2, sparse = FALSE))
(tr_function3 <- transition_matrix(Maze_function3, sparse = FALSE))

correct <- tr
expect_equal(tr, correct)
expect_equal(tr_dense, correct)
expect_equal(tr_sparse, correct)
expect_equal(tr_function2, correct)
expect_equal(tr_function3, correct)

# sparse get action matrix
tr <- transition_matrix(Maze_orig, 1, sparse = TRUE)
tr_dense <- transition_matrix(Maze_dense, 1,  sparse = TRUE)
tr_sparse <- transition_matrix(Maze_sparse, 1, sparse = TRUE)
tr_function2 <- transition_matrix(Maze_function2, 1, sparse = TRUE)
tr_function3 <- transition_matrix(Maze_function3, 1, sparse = TRUE)

correct <- tr
expect_equal(tr, correct)
expect_equal(tr_dense, correct)
expect_equal(tr_sparse, correct)
expect_equal(tr_function2, correct)
expect_equal(tr_function3, correct)

### action/row
(tr <- transition_matrix(Maze_orig, "up", 1))
(tr_dense <- transition_matrix(Maze_dense, "up", 1))
(tr_sparse <- transition_matrix(Maze_sparse, "up", 1))
(tr_function2 <- transition_matrix(Maze_function2, "up", 1))
(tr_function3 <- transition_matrix(Maze_function3, "up", 1))

correct <- tr
expect_equal(tr, correct)
expect_equal(tr_dense, correct)
expect_equal(tr_sparse, as(correct, "sparseVector"))
expect_equal(tr_function2, correct)
expect_equal(tr_function3, correct)


### action/col
(tr <- transition_matrix(Maze_orig, "up", , 2))
(tr_dense <- transition_matrix(Maze_dense, "up", , 2))
(tr_sparse <- transition_matrix(Maze_sparse, "up", , 2))
(tr_function2 <- transition_matrix(Maze_function2, "up", , 2))
(tr_function3 <- transition_matrix(Maze_function3, "up", , 2))

correct <- tr
expect_equal(tr, correct)
expect_equal(tr_dense, correct)
expect_equal(tr_sparse, as(correct, "sparseVector"))
expect_equal(tr_function2, correct)
expect_equal(tr_function3, correct)

### action/row/col
(tr <- transition_matrix(Maze_orig, "up", 1, 1))
(tr_dense <- transition_matrix(Maze_dense, "up", 1, 1))
(tr_sparse <- transition_matrix(Maze_sparse, "up", 1, 1))
(tr_function2 <- transition_matrix(Maze_function2, "up", 1, 1))
(tr_function3 <- transition_matrix(Maze_function3, "up", 1, 1))

correct <- tr
expect_equal(tr, correct)
expect_equal(tr_sparse, correct)
expect_equal(tr_dense, correct)
expect_equal(tr_function2, correct)
expect_equal(tr_function3, correct)

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


