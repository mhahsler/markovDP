# sample_sparse



p <- replicate(10, .1)

expect_equal(round_stochastic(p), p)
expect_equal(round_stochastic(p + 1e-9), p)

expect_true(sum1(p))
expect_true(sum1(round_stochastic(p + 1e-9)))
expect_true(sum1(round_stochastic(p + 1e-3)))
expect_true(sum1(p + 1e-9, digits = 7))
expect_false(sum1(p + 1e-9, digits = 10))
