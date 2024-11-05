# this is borrowed from nnet::which.is.max
which.max.random <- function(x, na.rm = FALSE) {
  mx <- x == max(x, na.rm = na.rm)
  mx[is.na(mx)] <- FALSE
  y <- seq_along(x)[mx]
  if (length(y) > 1L) {
    sample(y, 1L)
  } else {
    y
  }
}
