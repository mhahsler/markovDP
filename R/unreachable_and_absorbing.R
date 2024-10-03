#' Unreachable and Absorbing States
#'
#' Find unreachable and absorbing states using the transition model.
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
#' The function `absorbing_states()` checks if a state or a set of states are
#' absorbing (terminal states). A state is absorbing
#' if there is for all actions a probability of 1 for staying in the state.
#'
#' @name unreachable_and_absorbing
#' @aliases unreachable_and_absorbing
#' @family MDP
#'
#' @param model a [MDP] object.
#' @param state logical; check a single state. This can be much faster if
#'   the model contains a transition model implemented as a function.
#' @param horizon only states that can be reached within the
#'    horizon are reachable.
#' @param sparse logical; return a sparse logical vector?
#' @param use_precomputed logical; should precomputed values in the MDP be used?
#' @param progress logical; show a progress bar?
#' @param ... further arguments are passed on.
#'
#' @author Michael Hahsler
#' @examples
#' data(Maze)
#'
#' gridworld_matrix(Maze)
#' gridworld_matrix(Maze, what = "labels")
#'
#' # -1 and +1 are absorbing states
#' absorbing_states(Maze)
#' absorbing_states(Maze, sparse = FALSE)
#' absorbing_states(Maze, sparse = "states")
#'
#' # all states in the model are reachable
#' unreachable_states(Maze)
#' unreachable_states(Maze, sparse = FALSE)
#' unreachable_states(Maze, sparse = "states")
#' @importFrom Matrix colSums
#' @importFrom fastmap fastmap faststack
NULL

#' @rdname unreachable_and_absorbing
#' @returns  `unreachable_states()` returns a
#'  logical vector indicating the unreachable states.
#' @export
unreachable_states <- function(model,
                               horizon = Inf,
                               sparse = NULL,
                               use_precomputed = TRUE,
                               progress = TRUE,
                               ...) {
  UseMethod("unreachable_states")
}

#' @export
unreachable_states.MDP <- function(model,
                                   horizon = Inf,
                                   sparse = NULL,
                                   use_precomputed = TRUE,
                                   progress = TRUE,
                                   ...) {
  if (use_precomputed &&
      !is.null(model$unreachable_states) && !is.finite(horizon)) {
    return(.translate_logical(model$unreachable_states, model$states, sparse = sparse))
  }
  
  if (is.null(sparse))
    sparse <- TRUE
  
  .translate_logical(as(
    !.reachable_states(model, horizon = horizon, progress = progress),
    "sparseVector"
  ),
  model$states,
  sparse = sparse)
}


.reachable_states <- function(model,
                              horizon = Inf,
                              progress = TRUE) {
  reached <- fastmap()
  frontier <- faststack()
  
  for (start_state in start_vector(model, sparse = "states")) {
    frontier$push(start_state)
    # key: state label; value: depth
    reached$set(start_state, 0L)
  }
  
  if (progress) {
    pb <- my_progress_bar(N = length(model$states), name = "unreachable_states")
    pb$tick(0)
  }
  
  while (frontier$size() > 0) {
    state <- frontier$pop()
    
    if (progress) {
      pb$tick()
    }
    
    # available_actions is slow!
    #for (action in available_actions(model, state)){
    for (action in model$actions) {
      next_states <- transition_matrix(model, action, state, sparse = "states")
      
      for (next_state in next_states) {
        if (reached$has(next_state))
          next()
        
        depth <- reached$get(state)
        if (depth >= horizon)
          next()
        
        frontier$push(next_state)
        reached$set(next_state, depth + 1L)
      }
    }
  }
  
  if (progress)
    pb$terminate()
  
  model$states %in% names(reached$as_list())
}



#' @rdname unreachable_and_absorbing
#' @returns  `absorbing_states()` returns a logical vector indicating
#'    if the states are absorbing (terminal).
#' @export
absorbing_states <- function(model,
                             state = NULL,
                             sparse = NULL,
                             use_precomputed = TRUE,
                             ...) {
  UseMethod("absorbing_states")
}

#' @export
absorbing_states.MDP <- function(model,
                                 state = NULL,
                                 sparse = NULL,
                                 use_precomputed = TRUE,
                                 ...) {
  # do it faster to check a single state
  if (!is.null(state) &&
      length(state) == 1L) {
    if (use_precomputed && !is.null(model$absorbing_states)) {
      if (is.character(model$absorbing_states)) {
        if (!is.character(state))
          state <- model$states[state]
        return(state %in% model$absorbing_states)
      }
      # must be a logical vector (for sparse we need a index)
      if (is.character(state))
        state <- match(state, model$states)
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
      return(.translate_logical(model$absorbing_states, model$states, sparse))
    } else {
      absorbing <- rowSums(sapply(transition_matrix(model, sparse = NULL), diag)) == length(model$actions)
      if (is.null(sparse))
        sparse <- "states"
      return(.translate_logical(absorbing, model$states, sparse))
    }
  }
  
  # some states
  if (is.character(state))
    state <- match(state, model$states)
  
  if (use_precomputed && !is.null(model$absorbing_states)) {
    absorbing <- .translate_logical(model$absorbing_states, model$states, sparse = TRUE)
  } else {
    absorbing <- rowSums(sapply(transition_matrix(model, sparse = NULL), diag)) == length(model$actions)
  }
  
  absorbing <- absorbing[state]
  .translate_logical(absorbing, model$states[state], sparse)
}

#' @rdname unreachable_and_absorbing
#' @returns `remove_unreachable_states()` returns a model with all
#'  unreachable states removed.
#' @export
remove_unreachable_states <- function(model, ...) {
  keep <- !unreachable_states(model, sparse = FALSE)
  if (all(keep)) {
    return(model)
  }
  
  keep_states <- function(field, states) {
    if (is.data.frame(field)) {
      keep_names <- names(which(states))
      field <-
        field[field$start.state %in% c(NA, keep_names) &
                field$end.state %in% c(NA, keep_names), , drop = FALSE]
      field$start.state <-
        factor(as.character(field$start.state), levels = keep_names)
      field$end.state <-
        factor(as.character(field$end.state), levels = keep_names)
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
              m <- m[states, states, drop = FALSE]
            }
            m
          }
        )
    }
    field
  }
  
  
  # fix start state
  if (is.numeric(model$start) || is(model$start, "sparseVector")) {
    if (length(model$start) == length(model$states)) {
      ### prob vector
      model$start <- model$start[keep]
      if (!sum1(model$start)) {
        stop(
          "Probabilities for reachable states do not sum up to one! An unreachable state had a non-zero probability."
        )
      }
    } else {
      ### state ids... we translate to state names
      model$start <- model$states[model$start]
    }
  }
  if (is.character(model$start)) {
    if (model$start == "uniform") {
      # do nothing
    } else {
      model$start <- intersect(model$start, model$states[keep])
    }
    if (length(model$start) == 0L) {
      stop("Start state is not reachable.")
    }
  }
  
  model$states <- model$states[keep]
  model$transition_prob <- keep_states(model$transition_prob, keep)
  model$reward <- keep_states(model$reward, keep)
  
  # update reachable and unreachable
  if (!is.null(model$absorbing_states))
    model$absorbing_states <- model$absorbing_states[keep]
  if (!is.null(model$unreachable_states))
    model$unreachable_states <- model$unreachable_states[keep]
  
  # just check
  check_and_fix_MDP(model)
  model
}
