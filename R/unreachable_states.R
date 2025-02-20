#' Unreachable States
#'
#' Find or removes unreachable states using the transition model.
#'
#' The function `unreachable_states()` checks if states
#' cannot be reached from any other state. It performs a depth-first
#' search which can be slow.
#' The search breaks cycles to avoid an infinite loop.
#' The search depth can be restricted using
#' `horizon`.
#'
#' The function `remove_unreachable_states()` simplifies a model by
#' removing unreachable states from the model description.
#' 
#' @family MDP
#'
#' @param model a [MDP] object.
#' @param horizon only states that can be reached within the
#'    horizon are reachable.
#' @param sparse logical; return a sparse logical vector?
#' @param progress logical; show a progress bar?
#' @param ... further arguments are passed on.
#' 
#' @returns  `unreachable_states()` returns a
#'  logical vector indicating the unreachable states.
#'
#' @author Michael Hahsler
#' @examples
#' # create a Maze with an unreachable state
#' 
#' maze_unreach <- gw_read_maze(
#'     textConnection(c("XXXXXX", 
#'                      "XS X X",
#'                      "X  XXX",
#'                      "X   GX",
#'                      "XXXXXX")))
#' gw_plot(maze_unreach)
#'
#' unreachable_states(maze_unreach)
#' unreachable_states(maze_unreach, sparse = FALSE)
#' 
#' maze <- remove_unreachable_states(maze_unreach)
#' unreachable_states(maze)
#' gw_plot(maze)
#' @importFrom Matrix colSums
#' @export
unreachable_states <- function(model,
                                   horizon = Inf,
                                   sparse = "states",
                                   progress = TRUE,
                                   ...) {
  sparse <- sparse %||% TRUE
  if(is.null(S))
    stop("Model needs a defined state space.")
  
  states <- reachable_states(model, horizon = horizon, progress = progress, ...)
  unreachable <- setdiff(S(model), states)
  .translate_logical(unreachable,
                     S(model),
                     sparse = sparse)
}

#' @rdname unreachable_states
#' @returns `remove_unreachable_states()` returns a model with all
#'  unreachable states removed.
#' @export
remove_unreachable_states <- function(model, ...) {
  # Use dense because we need a named vector!
  unreachable <-  unreachable_states(model, sparse = "index")
  
  if (length(unreachable) == 0L) {
    model$unreachable <- character(0)
    return(model)
  }
  
  # remove states (ids) from a list of matrices or a data.frame
  remove_states <- function(field, states) {
    if (is.data.frame(field)) {
      field <-
        field[!(field$start.state %in% S(model)[states] |
                field$end.state %in% S(model)[states]), , drop = FALSE]
    } else if (is.function(field)) {
      # do nothing
    } else {
      ### a list of actions
      field <-
        lapply(
          field,
          FUN = function(m) {
            if (!is.character(m)) {
              ### strings like "uniform"
              m <- m[-states, -states, drop = FALSE]
            }
            m
          }
        )
    }
    field
  }
  
  
  # fix start state
  if (is.numeric(model$start) || is(model$start, "sparseVector")) {
    if (length(model$start) == length(S(model))) {
      ### prob vector
      model$start <- model$start[-unreachable]
      if (!sum1(model$start)) {
        stop(
          "Probabilities for reachable states do not sum up to one! An unreachable state had a non-zero probability."
        )
      }
    } else {
      ### state ids... we translate to state names and use code below!
      model$start <- setdiff(model$start, unreachable)
    }
  }
  if (is.character(model$start)) {
    if (model$start == "uniform") {
      # do nothing
    } else {
      model$start <- setdiff(model$start, S(model)[unreachable])
    }
    if (length(model$start) == 0L) {
      stop("Start state is not reachable.")
    }
  }
  
  # update reachable and unreachable
  if (!is.null(model$absorbing_states))
    model$absorbing_states <- setdiff(absorbing_states(model, sparse = "states"), S(model)[unreachable])
  
  model$states <- model$states[-unreachable]
  model$transition_prob <- remove_states(model$transition_prob, unreachable)
  model$reward <- remove_states(model$reward, unreachable)
  
  
  # just check
  check_and_fix_MDP(model)
  model
}
