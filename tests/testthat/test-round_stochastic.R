# round_stochastic a stochastic vector

# round a probability vector
p <- replicate(10, .1)

expect_equal(round_stochastic(p), p)
expect_equal(round_stochastic(p + 1e-9), p)
expect_equal(sum(round_stochastic(p + 1e-9)), 1)

# round a stochastic matrix
m <- matrix(runif(15), ncol = 3)
m <- sweep(m, 1, rowSums(m), "/")

expect_true(all(rowSums(round_stochastic(m, digits = 2)) == 1))
expect_true(all(rowSums(round_stochastic(m, digits = 1)) == 1))
expect_true(all(rowSums(round_stochastic(m, digits = 0)) == 1))

# sum1
expect_true(sum1(p))
expect_true(sum1(round_stochastic(p + 1e-9)))
expect_true(sum1(round_stochastic(p + 1e-3)))
expect_true(sum1(p + 1e-9, digits = 7))
expect_false(sum1(p + 1e-9, digits = 10))


