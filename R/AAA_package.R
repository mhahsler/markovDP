#' @keywords internal
#'
#' @section Key functions:
#' - Problem specification: [MDP]
#' - Solvers: [solve_MDP()]
#'
#' @import Rcpp
#' @importFrom Matrix spMatrix sparseVector crossprod coerce Math Math2 t diag rowSums cBind rBind nnzero
#' @importFrom methods as is new
#' @importFrom utils head tail read.table type.convert
#' @importFrom foreach foreach times %do% %dopar% getDoParWorkers
#' @useDynLib markovDP, .registration=TRUE
"_PACKAGE"

.onAttach <- function(libname, pkgname) {
  ### silence warning for no backend
  if (!foreach::getDoParRegistered()) {
    foreach::registerDoSEQ()
  }
}
