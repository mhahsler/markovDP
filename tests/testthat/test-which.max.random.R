# which.max.random

p <- replicate(4, 0.25)

# which always picks the first
expect_equal(length(table(replicate(10, which.max(p)))), 1L)

# should pick each
expect_equal(length(table(replicate(1000, which.max.random(p)))), 4L)
