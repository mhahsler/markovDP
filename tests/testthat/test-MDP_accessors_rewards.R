# Maze_orig contains a data.frame
# TODO: Add tests for functions

P <- transition_matrix(Maze_dense)

# dense
(r_dense <- reward_matrix(Maze_dense, sparse = FALSE))
(r_sparse <- reward_matrix(Maze_sparse, sparse = FALSE))
expect_true(inherits(r_dense[[1]], "matrix"))
expect_equal(r_dense, r_sparse)

# only compare for P>0
(r <- reward_matrix(Maze_dense, sparse = FALSE))
r <- mapply("*", lapply(P, ">", 0), r, SIMPLIFY = FALSE)
expect_equal(r, r_dense)

# sparse
(r_dense <- reward_matrix(Maze_dense, 1, sparse = TRUE))
(r_sparse <- reward_matrix(Maze_sparse, 1, sparse = TRUE))
expect_s4_class(r_dense, "dgCMatrix")
expect_equal(r_dense, r_sparse)

(r <- reward_matrix(Maze_dense, 1, sparse = TRUE))
non_zeroes <-  which(P[[1]] > 0, arr.ind = TRUE)
expect_s4_class(r, "dgCMatrix")
expect_equal(r[non_zeroes] , r_dense[non_zeroes])


### action/row (we use sparse = NULL)
(r_dense <- reward_matrix(Maze_dense, "up", 1, sparse = FALSE))
(r_sparse <- reward_matrix(Maze_sparse, "up", 1, sparse = FALSE))

expect_true(is.vector(r_dense))
expect_equal(r_dense, r_sparse)

(r <- reward_matrix(Maze_dense, "up", 1, sparse = FALSE))   ## should be df -> sparse
expect_true(is.vector(r))
non_zeroes <-  which(P[["up"]][1, ] > 0)
expect_equal(r_dense[non_zeroes], r[non_zeroes])

### action/row/col
(r <- reward_matrix(Maze_dense, "up", 1, 1))
(r_dense <- reward_matrix(Maze_dense, "up", 1, 1))
(r_sparse <- reward_matrix(Maze_sparse, "up", 1, 1))

expect_length(r, 1)
expect_equal(r, r_sparse)
expect_equal(r, r_dense)


# No simplify (when action is missing)

models <- c(models_matrix, list(Maze_orig, Maze_reward_function))

res <- models[[1]]$reward
for (m in models)
  expect_equal(reward_matrix(m, sparse = FALSE), res)

res <- lapply(models[[1]]$reward, "[", 1, 1)
for (m in models)
  expect_equal(reward_matrix(m, NULL, 1, 1), res)

res <- lapply(models[[1]]$reward, "[", 1:2, 1, drop = FALSE)
for (m in models)
  expect_equal(reward_matrix(m, NULL, 1:2, 1, sparse = FALSE), res)

res <- lapply(models[[1]]$reward, "[", 1:2, , drop = FALSE)
for (m in models)
  expect_equal(reward_matrix(m, NULL, 1:2, NULL, sparse = FALSE), res)

res <- lapply(models[[1]]$reward, "[", 1, 1:2, drop = FALSE)
for (m in models)
  expect_equal(reward_matrix(m, NULL, 1, 1:2, sparse = FALSE), res)

res <- lapply(models[[1]]$reward, "[", , 1:2, drop = FALSE)
for (m in models)
  expect_equal(reward_matrix(m, NULL, NULL, 1:2, sparse = FALSE), res)

res <- lapply(models[[1]]$reward, "[", 1:2, 1:2, drop = FALSE)
for (m in models)
  expect_equal(reward_matrix(m, NULL, 1:2, 1:2, sparse = FALSE), res)


# With simplify

# simplify does not fo anything if matrices are returned!
res <- models[[1]]$reward
for (m in models)
  expect_equal(reward_matrix(m, simplify = TRUE, sparse = FALSE), res)

res <- sapply(models[[1]]$reward, "[", 1, 1)
expect_true(is.vector(res))
for (m in models)
  expect_identical(reward_matrix(m, NULL, 1, 1, simplify = TRUE), res)

res <- sapply(models[[1]]$reward, "[", 1:2, 1)
expect_true(inherits(res, "matrix"))
for (m in models)
  expect_identical(reward_matrix(m, NULL, 1:2, 1, simplify = TRUE, sparse = FALSE), res)

res <- t(sapply(models[[1]]$reward, "[", 1, 1:2))
expect_true(inherits(res, "matrix"))
for (m in models)
  expect_identical(reward_matrix(m, NULL, 1, 1:2, simplify = TRUE, sparse = FALSE), res)

# simplify does not fo anything if matrices are returned!
res <- reward_matrix(models[[1]], NULL, 1:2, 1:2, simplify = TRUE, sparse = FALSE)
expect_true(inherits(res, "list"))
for (m in models)
  expect_identical(reward_matrix(m, NULL, 1:2, 1:2, simplify = TRUE, sparse = FALSE),
                   res)


# If trans is a function, reward is returned even if P = 0!

models <- c(models_trans_function, list(Maze_reward_trans_function))

rew <- reward_matrix(models[[1]], sparse = FALSE)

res <- rew
for (m in models)
  expect_equal(reward_matrix(m, sparse = FALSE), res)

res <- lapply(rew, "[", 1, 1)
for (m in models)
  expect_equal(reward_matrix(m, NULL, 1, 1), res)

res <- lapply(rew, "[", 1:2, 1, drop = FALSE)
for (m in models)
  expect_equal(reward_matrix(m, NULL, 1:2, 1, sparse = FALSE), res)

res <- lapply(rew, "[", 1:2, , drop = FALSE)
for (m in models)
  expect_equal(reward_matrix(m, NULL, 1:2, NULL, sparse = FALSE), res)

res <- lapply(rew, "[", 1, 1:2, drop = FALSE)
for (m in models)
  expect_equal(reward_matrix(m, NULL, 1, 1:2, sparse = FALSE), res)

res <- lapply(rew, "[", , 1:2, drop = FALSE)
for (m in models)
  expect_equal(reward_matrix(m, NULL, NULL, 1:2, sparse = FALSE), res)

res <- lapply(rew, "[", 1:2, 1:2, drop = FALSE)
for (m in models)
  expect_equal(reward_matrix(m, NULL, 1:2, 1:2, sparse = FALSE), res)


