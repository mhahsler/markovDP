# Maze_orig contains a data.frame

# dense
(tr <- reward_matrix(Maze_orig, sparse = FALSE))
(tr_dense <- reward_matrix(Maze_dense, sparse = FALSE))
(tr_sparse <- reward_matrix(Maze_sparse, sparse = FALSE))

correct <- tr

#expect_equal(tr, correct)
expect_equal(tr_dense, correct)

# sparse representation has 0 where P == 0
P <- transition_matrix(Maze_dense)
correct_sparse <- sapply(A(Maze_dense), 
                         FUN = function(a) { ca <- correct[[a]]; ca[P[[a]] == 0] <- 0; ca },
                        simplify = FALSE)
expect_equal(tr_sparse, correct_sparse)

# sparse
(tr <- reward_matrix(Maze_orig, 1, sparse = TRUE))
(tr_dense <- reward_matrix(Maze_dense, 1,  sparse = TRUE))
(tr_sparse <- reward_matrix(Maze_sparse, 1, sparse = TRUE))

correct <- tr
#expect_equal(tr, correct)
expect_equal(tr_dense, correct)
expect_equal(tr_sparse, as(correct_sparse[[1]], "dgCMatrix"))

### action/row (we use sparse = NULL)
(tr <- reward_matrix(Maze_orig, "up", 1, sparse = FALSE))   ## should be df -> sparse
(tr_dense <- reward_matrix(Maze_dense, "up", 1, sparse = FALSE))
(tr_sparse <- reward_matrix(Maze_sparse, "up", 1, sparse = FALSE))

correct <- tr
#expect_equal(tr, correct)
expect_equal(tr_dense, correct)
expect_equal(tr_sparse, correct_sparse[[1]][1, ])

### action/row/col
(tr <- reward_matrix(Maze_orig, "up", 1, 1))
(tr_dense <- reward_matrix(Maze_dense, "up", 1, 1))
(tr_sparse <- reward_matrix(Maze_sparse, "up", 1, 1))

correct <- tr
#expect_equal(tr, correct)
expect_equal(tr_sparse, correct)
expect_equal(tr_dense, correct)


# simplify (when action is missing)
for (m in models) {
  (tr <- reward_matrix(m, NULL, 1, 1))
  expect_true(inherits(tr, "list"))
  (tr <- reward_matrix(m, NULL, 1:2, 1))
  expect_true(inherits(tr, "list"))
  (tr <- reward_matrix(m, NULL, 1, 1:2))
  expect_true(inherits(tr, "list"))
  (tr <- reward_matrix(m, NULL, 1, 1:2))
  expect_true(inherits(tr, "list"))
  (tr <- reward_matrix(m, NULL, 1:2, 1:2))
  expect_true(inherits(tr, "list"))
}

res <- sapply(models[[1]]$reward, "[", 1, 1)
# reward_matrix(models[[1]], NULL, 1, 1, simplify = TRUE)
expect_true(inherits(res, "numeric"))
for (m in models)
  expect_identical(reward_matrix(m, NULL, 1, 1, simplify = TRUE), res)

res <- sapply(models[[1]]$reward, "[", 1:2, 1)
#  reward_matrix(models[[1]], NULL, 1:2, 1, simplify = TRUE, sparse = FALSE)
expect_true(inherits(res, "matrix"))
for (m in models)
  expect_identical(reward_matrix(m, NULL, 1:2, 1, simplify = TRUE, sparse = FALSE), res)

res <- t(sapply(models[[1]]$reward, "[", 1, 1:2)) 
#reward_matrix(models[[1]], NULL, 1, 1:2, simplify = TRUE, sparse = FALSE)
expect_true(inherits(res, "matrix"))
for (m in models[-2])
  expect_identical(reward_matrix(m, NULL, 1, 1:2, simplify = TRUE, sparse = FALSE), res)

res <- reward_matrix(models[[1]], NULL, 1:2, 1:2, simplify = TRUE, sparse = FALSE)
expect_true(inherits(res, "list"))
for (m in models[-2])
  expect_identical(reward_matrix(m, NULL, 1:2, 1:2, simplify = TRUE, 
                                     sparse = FALSE), res)



