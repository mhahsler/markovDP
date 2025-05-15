near_equal <- function(x, y, tol = .Machine$double.eps) 
  abs(x - y) < tol 