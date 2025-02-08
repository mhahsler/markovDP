#' Sample a Start State
#'
#' Samples a start state which can be used for simulation.
#' 
#' @family MDP
#' @family MDPTF
#'
#' @param model an MDP model.
#'
#' @returns a start state.
#'
#' @author Michael Hahsler
#' @examples
#' data(Maze)
#'
#' start(Maze)
#' @export
start <- function(model) {
  UseMethod("start")
}

#' @export
start.MDP <- function(model) {
  if (length(model$start) == 1L)
    s <- model$start
  else
    s <- sample.int(length(S(model)), 1L, prob = start_vector(model, sparse = FALSE))
  normalize_state(s, model)
}
  
#' @export
start.MDPTF <- function(model) {
  model$start
}