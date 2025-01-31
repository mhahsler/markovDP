#' Absorbing States
#'
#' Find absorbing states using the transition model.
#'
#' The function `absorbing_states()` checks if a state or a set of states are
#' absorbing (terminal states). A state is absorbing
#' if there is for all actions a probability of 1 for staying in the state.
#'
#' @family MDP
#' @family MDPE
#'
#' @param model a [MDP] object.
#' @param state a single state to check. This can be much faster if
#'   the model contains a transition model implemented as a function.
#'   `NULL` means all states are checked.
#' @param sparse logical; if return a sparse logical vector?
#' @param use_precomputed logical; should precomputed values in the MDP be used?
#' @param ... further arguments are passed on.
#' 
#' @returns  `absorbing_states()` returns a logical vector indicating
#'    if the states are absorbing (terminal).
#'
#' @author Michael Hahsler
#' @examples
#' data(Maze)
#'
#' gw_matrix(Maze)
#' gw_matrix(Maze, what = "labels")
#' gw_matrix(Maze, what = "absorbing")
#'
#' # -1 and +1 are absorbing states
#' absorbing_states(Maze)
#' absorbing_states(Maze, sparse = FALSE)
#' absorbing_states(Maze, sparse = "states")
#' 
#' # check individual states
#' absorbing_states(Maze, "s(1,1)")
#' absorbing_states(Maze, "s(1,4)")
#' @importFrom Matrix colSums
#' @export
absorbing_states <- function(model,
                             state = NULL,
                             sparse = "states",
                             use_precomputed = TRUE,
                             ...) {
  UseMethod("absorbing_states")
}

#' @export
absorbing_states.MDP <- function(model,
                                 state = NULL,
                                 sparse = "states",
                                 use_precomputed = TRUE,
                                 ...) {
  # do it faster to check a single state
  if (!is.null(state) &&
      length(state) == 1L) {
    if (use_precomputed && !is.null(model$absorbing_states)) {
      if (is.character(model$absorbing_states)) {
        if (!is.character(state))
          state <- S(model)[state]
        return(state %in% model$absorbing_states)
      }
      # must be a logical vector (for sparse we need a index)
      if (is.character(state))
        state <- match(state, S(model))
      return(as(model$absorbing_states[state], "vector"))
    } else
      return(all(
        transition_matrix(
          model,
          start.state = state,
          end.state = state,
          simplify = TRUE
        ) == 1
      ))
  }
  
  # all states
  if (is.null(state)) {
    if (use_precomputed && !is.null(model$absorbing_states)) {
      return(.translate_logical(model$absorbing_states, S(model), sparse))
    } else {
      absorbing <- rowSums(sapply(transition_matrix(model, sparse = NULL), diag)) == length(A(model))
      if (is.null(sparse))
        sparse <- "states"
      return(.translate_logical(absorbing, S(model), sparse))
    }
  }
  
  # some states
  if (is.character(state))
    state <- match(state, S(model))
  
  if (use_precomputed && !is.null(model$absorbing_states)) {
    absorbing <- .translate_logical(model$absorbing_states, S(model), sparse = TRUE)
  } else {
    absorbing <- rowSums(sapply(transition_matrix(model, sparse = NULL), diag)) == length(A(model))
  }
  
  absorbing <- absorbing[state]
  .translate_logical(absorbing, S(model)[state], sparse)
}

#' @export
absorbing_states.MDPE <- function(model,
                                 state = NULL,
                                 sparse = "states",
                                 use_precomputed = TRUE,
                                 ...) {
  
  if(is.null(model$absorbing_states))
    stop("Model does not specify absorbing states!")
  
  if(is.null(state))
    return(model$absorbing_states)
  
  state <- drop(normalize_state_features(state, model))
  
  any(apply(model$absorbing_states, MARGIN = 1, FUN = function(x) all(x == state)))
}
