# sample that accepts dense or sparse probability vectors
sample_sparse <- function(x, size, replace = FALSE, prob) {
  if (is(prob, "sparseVector")) {
    x[prob@i][sample.int(length(prob@i), size, replace, prob = prob@x)]
  } else 
    x[sample.int(length(x), size, replace, prob = prob)]
}