.onLoad <- function(libname, pkgname){
  # Use faster dense matrices for small problems
  # Use about 100 MB with one double using 8 bytes 
  
  options("MDP_SPARSE_LIMIT" = 100e6 / 8)
}