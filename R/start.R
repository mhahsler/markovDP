#' Sample a Start State
#'
#' Samples a start state which can be used for simulation.
#' 
#' @family MDP
#' @family MDPTF
#'
#' @aliases start
#' @name start
#' @param x an MDP model.
#'
#' @returns a start state.
#'
#' @author Michael Hahsler
#' @examples
#' data(Maze)
#'
#' start(Maze)
#' @importFrom stats start
NULL

# generic is defined in stats
# @export
#start <- function(x, ...) {
#  UseMethod("start")
#}

#' @rdname start
#' @param x a model.
#' @param as a returned format. See parameter `as` in [`normalize_state()`].
#' @param ... additional parameters are ignored.
#' @export
start.MDP <- function(x, as = "factor", ...) {
  .nodots(...)
  model <- x
  if (length(model$start) == 1L)
    s <- model$start
  else
    s <- sample.int(length(S(model)), 1L, prob = start_vector(model, sparse = FALSE))
  normalize_state(s, model, as = as)
}
  
#' @rdname start
#' @export
start.MDPTF <- function(x, as = "features", ...) {
  .nodots(...)
  normalize_state(x$start, x, as = as)
}