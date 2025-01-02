verbose <- FALSE
# verbose <- TRUE

state_labels <- paste0("s", 1:10)

states_empty <- list(
  character(0), 
  integer(0),
  structure(c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
            names = state_labels),
  sparseVector(i = integer(0), x = logical(0), length = 10))

states_single <- list(
  "s2", 
  2,
  structure(c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
            names = state_labels),
  sparseVector(i = 2L, x = TRUE, length = 10))

states_single <- list( 
  c("s2", "s4"), 
  c(2, 4),
  structure(c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
            names = state_labels),
  sparseVector(i = c(2L, 4L), x = c(TRUE, TRUE), length = 10))

sparse_lvl <- list(
  "states",
  "index",
  FALSE,
  TRUE)


for (test_case in list(states_empty, states_single, states_single)) {
  for (s in test_case) {
    for (sl in seq_along(sparse_lvl)) { 
      
      if(verbose) {
        cat("--------------------\n")
        cat("input = "); print(s)
        cat("sparse = ", sparse_lvl[[sl]], "\n")
        cat("expected = "); print(test_case[[sl]]); cat("\n")
      }
      
      tr <- .translate_logical(s, state_labels, sparse = sparse_lvl[[sl]])
      
      if(verbose) {
        cat("output = "); print(tr); cat("\n")
      }
      
      expect_equal(tr, test_case[[sl]])
    }
  }  
}
