# dense
(sv_dense <- start_vector(Maze_orig, sparse = FALSE))
(sv_sparse <- start_vector(Maze_orig, sparse = TRUE))

expect_equal(sv_sparse, as(sv_dense, "sparseVector"))

# Keyword uniform
Maze_test <- Maze_orig
Maze_test$start <- "uniform"

(sv_dense <- start_vector(Maze_test, sparse = FALSE))
(sv_sparse <- start_vector(Maze_test, sparse = TRUE))
expect_equal(unname(sv_dense), 
             rep(1/length(Maze_test$states), length(Maze_test$states)))
expect_equal(sv_sparse, as(sv_dense, "sparseVector"))

# 3 states by id
Maze_test <- Maze_orig
Maze_test$start <- 2:4
correct <- numeric(length(Maze_test$states))
correct[2:4] <- 1/3

(sv_dense <- start_vector(Maze_test, sparse = FALSE))
(sv_sparse <- start_vector(Maze_test, sparse = TRUE))
expect_equal(unname(sv_dense), correct)
expect_equal(sv_sparse, as(correct, "sparseVector"))

# 3 states by name
Maze_test <- Maze_orig
Maze_test$start <- Maze_orig$states[2:4]
correct <- numeric(length(Maze_test$states))
correct[2:4] <- 1/3

(sv_dense <- start_vector(Maze_test, sparse = FALSE))
(sv_sparse <- start_vector(Maze_test, sparse = TRUE))
expect_equal(unname(sv_dense), correct)
expect_equal(sv_sparse, as(correct, "sparseVector"))

# 3 states removed by neg. id
Maze_test$start <- -1 * (2:4)
correct <- numeric(length(Maze_test$states))
correct[c(1, 5:11)] <- 1/8

(sv_dense <- start_vector(Maze_test, sparse = FALSE))
(sv_sparse <- start_vector(Maze_test, sparse = TRUE))
expect_equal(unname(sv_dense), correct)
expect_equal(sv_sparse, as(correct, "sparseVector"))

# 3 states removed by name
Maze_test$start <- c("-", Maze_orig$states[2:4])
correct <- numeric(length(Maze_test$states))
correct[c(1, 5:11)] <- 1/8

(sv_dense <- start_vector(Maze_test, sparse = FALSE))
(sv_sparse <- start_vector(Maze_test, sparse = TRUE))
expect_equal(unname(sv_dense), correct)
expect_equal(sv_sparse, as(correct, "sparseVector"))
