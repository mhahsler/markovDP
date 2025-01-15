# sample sparse

# create a sparse vector
non_zeros <- 10
#non_zeros <- 100
#non_zeros <- 500
p <- numeric(1000)
non_zero_idx <- sample.int(length(p), non_zeros)
p[non_zero_idx] <- runif(non_zeros)
p <- p/sum(p)

(p_sparse <- as(p, "sparseVector"))

# base sample does not need to be tested
#(samp <- sample.int(length(p), size = 100, replace = TRUE, prob = p))
#expect_true(all(samp %in% non_zero_idx))

(samp <- sample_sparse(seq_along(p), size = 100, replace = TRUE, prob = p_sparse))
expect_true(all(samp %in% non_zero_idx))


# sample.int is faster, but sparse sampling takes less memory!
#bench::mark(samp <- sample.int(length(p), size = 100, replace = TRUE, prob = p))
#bench::mark(samp <- sample_sparse(seq_along(p), size = 100, replace = TRUE, prob = p))
#bench::mark(samp <- sample_sparse(seq_along(p), size = 100, replace = TRUE, prob = p_sparse))

