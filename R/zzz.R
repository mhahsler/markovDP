.onLoad <- function(libname, pkgname){
  # Use faster dense matrices for small problems
  # Use about 1 GB with one double using 8 bytes 
  options("MDP_SPARSE_LIMIT" = 1e9 / 8)
  
  MatrixExtra::set_new_matrix_behavior()
}