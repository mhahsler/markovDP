# dense
(sv_dense <- start_vector(Maze_orig, sparse = FALSE))
(sv_sparse <- start_vector(Maze_orig, sparse = TRUE))

expect_equal(sv_sparse, as(sv_dense, "sparseVector"))

# Keyword uniform
model_test <- Maze_orig
model_test$start <- "uniform"

(sv_dense <- start_vector(model_test, sparse = FALSE))
(sv_sparse <- start_vector(model_test, sparse = TRUE))
expect_equal(unname(sv_dense), 
             rep(1/length(model_test$states), length(model_test$states)))
expect_equal(sv_sparse, as(sv_dense, "sparseVector"))

# 3 states by id
model_test <- Maze_orig
model_test$start <- 2:4
correct <- numeric(length(model_test$states))
correct[2:4] <- 1/3

(sv_dense <- start_vector(model_test, sparse = FALSE))
(sv_sparse <- start_vector(model_test, sparse = TRUE))
expect_equal(unname(sv_dense), correct)
expect_equal(sv_sparse, as(correct, "sparseVector"))

# 3 states by name
model_test <- Maze_orig
model_test$start <- Maze_orig$states[2:4]
correct <- numeric(length(model_test$states))
correct[2:4] <- 1/3

(sv_dense <- start_vector(model_test, sparse = FALSE))
(sv_sparse <- start_vector(model_test, sparse = TRUE))
expect_equal(unname(sv_dense), correct)
expect_equal(sv_sparse, as(correct, "sparseVector"))

# 3 states removed by neg. id
model_test$start <- -1 * (2:4)
correct <- numeric(length(model_test$states))
correct[c(1, 5:11)] <- 1/8

(sv_dense <- start_vector(model_test, sparse = FALSE))
(sv_sparse <- start_vector(model_test, sparse = TRUE))
expect_equal(unname(sv_dense), correct)
expect_equal(sv_sparse, as(correct, "sparseVector"))

# 3 states removed by name
model_test$start <- c("-", Maze_orig$states[2:4])
correct <- numeric(length(model_test$states))
correct[c(1, 5:11)] <- 1/8

(sv_dense <- start_vector(model_test, sparse = FALSE))
(sv_sparse <- start_vector(model_test, sparse = TRUE))
expect_equal(unname(sv_dense), correct)
expect_equal(sv_sparse, as(correct, "sparseVector"))
