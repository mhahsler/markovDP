# Maze_orig contains a data.frame

(R1_full <- reward_matrix(Maze_function2, 1, sparse = FALSE))
(R1_non_zero <- (Maze_dense$transition_prob[[1]] != 0) * R1_full)

# Translation rules:
# * if trans dense ->  dense (P == 0 missing).    Example: Maze_orig
# * if trans sparse ->  sparse (P == 0 missing)   
# * if trans not matrix (function/DF) -> sparse/dense (complete)
#                       Example: Maze_function2, Maze_reward_trans_function

### Full matrix
# sparse = NULL
(r_orig <- reward_matrix(Maze_orig, sparse = NULL))
expect_equal(r_orig[[1]], R1_non_zero)

(r_dense <- reward_matrix(Maze_dense, sparse = NULL))
expect_true(is.matrix(r_dense[[1]]))
expect_equal(r_dense[[1]], R1_non_zero)

(r_sparse <- reward_matrix(Maze_sparse, sparse = NULL))
expect_s4_class(r_sparse[[1]], "dgRMatrix")
expect_equal(as.matrix(r_sparse[[1]]), R1_non_zero)

# small problems -> dense
(r_func <- reward_matrix(Maze_function2, sparse = NULL))
(r_all_func <- reward_matrix(Maze_reward_trans_function, sparse = NULL))
expect_equal(r_func[[1]], R1_full)
expect_equal(r_func, r_all_func)

# sparse = FALSE
(r_orig <- reward_matrix(Maze_orig, sparse = FALSE))
(r_dense <- reward_matrix(Maze_dense, sparse = FALSE))
(r_sparse <- reward_matrix(Maze_sparse, sparse = FALSE))
expect_equal(r_sparse[[1]], R1_non_zero)
expect_equal(r_sparse, r_dense)
expect_equal(r_sparse, r_orig)

(r_func <- reward_matrix(Maze_function2, sparse = FALSE))
(r_all_func <- reward_matrix(Maze_reward_trans_function, sparse = FALSE))
expect_equal(r_func[[1]], R1_full)
expect_equal(r_func, r_all_func)

# sparse = TRUE
(r_orig <- reward_matrix(Maze_orig, sparse = TRUE))
(r_dense <- reward_matrix(Maze_dense, sparse = TRUE))
(r_sparse <- reward_matrix(Maze_sparse, sparse = TRUE))
expect_s4_class(r_sparse[[1]], "dgRMatrix")
expect_equal(as.matrix(r_sparse[[1]]), R1_non_zero)
expect_equal(r_sparse, r_dense)
expect_equal(r_sparse, r_orig)

(r_func <- reward_matrix(Maze_function2, sparse = TRUE))
(r_all_func <- reward_matrix(Maze_reward_trans_function, sparse = TRUE))
expect_s4_class(r_func[[1]], "dgRMatrix")
expect_equal(as.matrix(r_func[[1]]), R1_full)
expect_equal(r_func, r_all_func)

### just action is too simple to test

### action/row 
# sparse = NULL

(r_orig <- reward_matrix(Maze_orig, 1, 1, sparse = NULL))
expect_equal(r_orig, R1_non_zero[1,])

(r_dense <- reward_matrix(Maze_dense, 1, 1, sparse = NULL))
expect_true(is.vector(r_dense))
expect_equal(r_dense, R1_non_zero[1,])

(r_sparse <- reward_matrix(Maze_sparse, 1, 1, sparse = NULL))
expect_s4_class(r_sparse, "dsparseVector")
expect_equal(as.vector(r_sparse), unname(R1_non_zero[1,]))

# small problems -> dense
(r_func <- reward_matrix(Maze_function2, 1, 1, sparse = NULL))
(r_all_func <- reward_matrix(Maze_reward_trans_function, 1, 1, sparse = NULL))
expect_equal(r_func, R1_full[1, ])
expect_equal(r_func, r_all_func)

# sparse = FALSE
(r_orig <- reward_matrix(Maze_orig, 1, 1, sparse = FALSE))
expect_equal(r_orig, R1_non_zero[1,])

(r_dense <- reward_matrix(Maze_dense, 1, 1, sparse = FALSE))
(r_sparse <- reward_matrix(Maze_sparse, 1, 1, sparse = FALSE))
expect_equal(r_sparse, R1_non_zero[1, ])
expect_equal(r_sparse, r_dense)

(r_func <- reward_matrix(Maze_function2, 1, 1, sparse = FALSE))
(r_all_func <- reward_matrix(Maze_reward_trans_function, 1, 1, sparse = FALSE))
expect_equal(r_func, R1_full[1,])
expect_equal(r_func, r_all_func)

# sparse = TRUE
(r_orig <- reward_matrix(Maze_orig, 1, 1, sparse = TRUE))
(r_dense <- reward_matrix(Maze_dense, 1, 1, sparse = TRUE))
(r_sparse <- reward_matrix(Maze_sparse, 1, 1, sparse = TRUE))
expect_s4_class(r_sparse, "dsparseVector")
expect_equal(as.vector(r_sparse), unname(R1_non_zero[1, ]))
expect_equal(r_sparse, r_dense)
expect_equal(r_sparse, r_orig)

(r_func <- reward_matrix(Maze_function2, 1, 1, sparse = TRUE))
(r_all_func <- reward_matrix(Maze_reward_trans_function, 1, 1, sparse = TRUE))
expect_s4_class(r_func, "dsparseVector")
expect_equal(as.vector(r_func), unname(R1_full[1,]))
expect_equal(r_func, r_all_func)

### action/row/col
# all return the same
(r_orig <- reward_matrix(Maze_orig, 1, 1, 1))
(r_dense <- reward_matrix(Maze_dense, 1, 1, 1))
(r_sparse <- reward_matrix(Maze_sparse, 1, 1, 1))
(r_func <- reward_matrix(Maze_function2, 1, 1, 1))
(r_all_func <- reward_matrix(Maze_reward_trans_function, 1, 1, 1))

expect_length(r_dense, 1)
expect_equal(r_dense, r_orig)
expect_equal(r_dense, r_sparse)
expect_equal(r_dense, r_func)
expect_equal(r_dense, r_all_func)

# ranges
(r_orig <- reward_matrix(Maze_orig, 1, 1:2, 1:2))
(r_dense <- reward_matrix(Maze_dense, 1, 1:2, 1:2))
(r_sparse <- reward_matrix(Maze_sparse, 1, 1:2, 1:2))
(r_func <- reward_matrix(Maze_function2, 1, 1:2, 1:2))
(r_all_func <- reward_matrix(Maze_reward_trans_function, 1, 1:2, 1:2))
expect_equal(dim(r_dense), c(2L, 2L))
expect_equal(r_dense,  R1_non_zero[1:2,1:2])
expect_equal(r_dense, as.matrix(r_sparse))
expect_equal(r_dense, r_orig)

expect_equal(dim(r_func), c(2L, 2L))
expect_equal(r_func,  R1_full[1:2,1:2])
expect_equal(r_func, r_all_func)

(r_orig <- reward_matrix(Maze_orig, 1, 1:2, 1))
(r_dense <- reward_matrix(Maze_dense, 1, 1:2, 1))
(r_sparse <- reward_matrix(Maze_sparse, 1, 1:2, 1))
(r_func <- reward_matrix(Maze_function2, 1, 1:2, 1))
(r_all_func <- reward_matrix(Maze_reward_trans_function, 1, 1:2, 1))
expect_equal(dim(r_dense), c(2L, 1L))
expect_equal(r_dense,  R1_non_zero[1:2,1,drop = FALSE])
expect_equal(r_dense, as.matrix(r_sparse))
expect_equal(r_dense, r_orig)

expect_equal(dim(r_func), c(2L, 1L))
expect_equal(r_func,  R1_non_zero[1:2,1, drop = FALSE])
expect_equal(r_func, r_all_func)

(r_orig <- reward_matrix(Maze_orig, 1, 1, 1:2))
(r_dense <- reward_matrix(Maze_dense, 1, 1, 1:2))
(r_sparse <- reward_matrix(Maze_sparse, 1, 1, 1:2))
(r_func <- reward_matrix(Maze_function2, 1, 1, 1:2))
(r_all_func <- reward_matrix(Maze_reward_trans_function, 1, 1, 1:2))
expect_equal(dim(r_dense), c(1L, 2L))
expect_equal(r_dense,  R1_non_zero[1,1:2, drop = FALSE])
expect_equal(r_dense, as.matrix(r_sparse))
expect_equal(r_dense, r_orig)

expect_equal(dim(r_func), c(1L, 2L))
expect_equal(r_func,  R1_full[1,1:2, drop = FALSE])
expect_equal(r_func, r_all_func)

# action is missing
(r_orig <- reward_matrix(Maze_orig, NULL, 1, 1, ))
(r_dense <- reward_matrix(Maze_dense, NULL, 1, 1, ))
(r_sparse <- reward_matrix(Maze_sparse, NULL, 1, 1, ))
(r_func <- reward_matrix(Maze_function2, NULL, 1, 1))
(r_all_func <- reward_matrix(Maze_reward_trans_function, NULL, 1, 1))

expect_true(is.list(r_dense))
expect_equal(r_dense, r_orig)
expect_equal(r_dense, r_sparse)
expect_equal(r_dense, r_func)
expect_equal(r_dense, r_all_func)

(r_orig <- reward_matrix(Maze_orig, NULL, 1, 1, simplify = TRUE))
(r_dense <- reward_matrix(Maze_dense, NULL, 1, 1, simplify = TRUE))
(r_sparse <- reward_matrix(Maze_sparse, NULL, 1, 1, simplify = TRUE))
(r_func <- reward_matrix(Maze_function2, NULL, 1, 1, simplify = TRUE))
(r_all_func <- reward_matrix(Maze_reward_trans_function, NULL, 1, 1, simplify = TRUE))

expect_true(is.vector(r_dense))
expect_equal(r_dense, r_orig)
expect_equal(r_dense, r_sparse)
expect_equal(r_dense, r_func)
expect_equal(r_dense, r_all_func)

