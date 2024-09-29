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
#' absorbing_states(Maze, sparse = TRUE)
#' which(absorbing_states(Maze))
#'
#' # all states in the model are reachable
#' unreachable_states(Maze)
#' unreachable_states(Maze, sparse = TRUE)
#' which(unreachable_states(Maze))
#' @importFrom Matrix colSums
NULL

#' @rdname unreachable_and_absorbing
#' @returns  `unreachable_states()` returns a
#'  logical vector indicating the unreachable states.
#' @export
unreachable_states <- function(model,
                               horizon = Inf,
                               sparse = FALSE,
                               progress = TRUE,
                               ...) {
  UseMethod("unreachable_states")
}

#' @export
unreachable_states.MDP <- function(model,
                                   horizon = Inf,
                                   sparse = FALSE,
                                   progress = TRUE,
                                   ...) {
  if (!is.null(model$unreachable_states) && !is.finite(horizon)) {
    return(
      .sparsify_vector(
        model$unreachable_states,
        sparse = sparse,
        names = model$states
      )
    )
  }

 .sparsify_vector(!.reachable_states(model, 
                                     horizon = horizon, 
                                     progress = progress), 
                  sparse = sparse, 
                  names = model$states)
}


.reachable_states <- function(model, horizon = Inf, progress = TRUE) {
  states <- new.env(hash = TRUE)
  
  if (progress) {
    pb <- my_progress_spinner(name = "unreachable_states")
    pb$tick(0)
  }
  
  move <- function(state, depth = 0) {
    if (depth > horizon)
      return()
    
    for (action in available_actions(model, state)){
      next_states <- transition_matrix(model, action, state, sparse = "states")  
      
      for (next_state in next_states) {
        if (exists(next_state, envir = states, inherits = FALSE))
          next()
        
        if (progress)
          pb$tick()
        
        assign(next_state, TRUE, envir = states)
        move(next_state, depth + 1L)
      }
    }
  }
  
  for (start_state in start_vector(model, sparse = "states")) {
    assign(start_state, TRUE, envir = states)
    move(start_state)
  }
  
  Vectorize(exists, vectorize.args = "x")(model$states, 
                                          envir = states, 
                                          inherits = FALSE)
}



#' @rdname unreachable_and_absorbing
#' @returns  `absorbing_states()` returns a logical vector indicating
#'    if the states are absorbing (terminal).
#' @export
absorbing_states <- function(model,
                             state = NULL,
                             sparse = FALSE,
                             ...) {
  UseMethod("absorbing_states")
}

#' @export
absorbing_states.MDP <- function(model,
                                 state = NULL,
                                 sparse = FALSE,
                                 ...) {
  
  # do it faster to check a single state
  if (!is.null(state) &&
      length(state) == 1L) { 
    if (!is.null(model$absorbing_states)) {
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
      
  if (!is.null(model$absorbing_states)) {
    absorbing <- .sparsify_vector(model$absorbing_states, sparse = sparse, names = model$states)
  } else {
    absorbing <- rowSums(sapply(transition_matrix(model, sparse = NULL), diag)) == length(model$actions)
    absorbing <- .sparsify_vector(absorbing, sparse, names = model$states)
  }
  
  if (!is.null(state)) {
    if (is.character(state))
      state <- match(state, model$states)
    absorbing <- absorbing[state]
  }
  
  absorbing
}

#' @rdname unreachable_and_absorbing
#' @returns `remove_unreachable_states()` returns a model with all
#'  unreachable states removed.
#' @export
remove_unreachable_states <- function(model, ...) {
  reachable <- !unreachable_states(model, ...)
  if (all(reachable)) {
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
      model$start <- model$start[reachable]
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
      model$start <- intersect(model$start, model$states[reachable])
    }
    if (length(model$start) == 0L) {
      stop("Start state is not reachable.")
    }
  }
  
  model$states <- model$states[reachable]
  model$transition_prob <- keep_states(model$transition_prob, reachable)
  model$reward <- keep_states(model$reward, reachable)
  if (!is.null(model$observations)) {
    model$observations <- keep_states(model$observations, reachable)
  }
  
  # update reachable and unreachable
  if (!is.null(model$absorbing_states))
    model$absorbing_states <- model$absorbing_states[reachable]
  if (!is.null(model$reachable_states))
    model$reachable_states <- model$reachable_states[reachable]
  
  # just check
  check_and_fix_MDP(model)
  model
}
