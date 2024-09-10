# Maze_orig contains a data.frame

# dense
(tr <- reward_matrix(Maze_orig, sparse = FALSE))
(tr_dense <- reward_matrix(Maze_dense, sparse = FALSE))
(tr_sparse <- reward_matrix(Maze_sparse, sparse = FALSE))

correct <- tr

#expect_equal(tr, correct)
expect_equal(tr_dense, correct)
expect_equal(tr_sparse, correct)

# sparse
(tr <- reward_matrix(Maze_orig, 1, sparse = TRUE))
(tr_dense <- reward_matrix(Maze_dense, 1,  sparse = TRUE))
(tr_sparse <- reward_matrix(Maze_sparse, 1, sparse = TRUE))

correct <- tr
#expect_equal(tr, correct)
expect_equal(tr_dense, correct)
expect_equal(tr_sparse, correct)

### action/row (we use sparse = NULL)
(tr <- reward_matrix(Maze_orig, "up", 1))   ## should be df -> sparse
(tr_dense <- reward_matrix(Maze_dense, "up", 1))
(tr_sparse <- reward_matrix(Maze_sparse, "up", 1))

correct <- as(tr, "sparseVector") ## make sure it is sparse
expect_equal(tr, correct)
expect_equal(as(tr_dense, "sparseVector") , correct)
expect_equal(tr_sparse, correct)

### action/row/col
(tr <- reward_matrix(Maze_orig, "up", 1, 1))
(tr_dense <- reward_matrix(Maze_dense, "up", 1, 1))
(tr_sparse <- reward_matrix(Maze_sparse, "up", 1, 1))

correct <- tr
#expect_equal(tr, correct)
expect_equal(tr_sparse, correct)
expect_equal(tr_dense, correct)
