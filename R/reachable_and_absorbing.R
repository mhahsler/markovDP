#' Unreachable and Absorbing States
#'
#' Find unreachable and absorbing states using the transition model.
#'
#' The function `unreachable_states()` checks if states
#' cannot be reached from any other state.
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
#' @param sparse logical; return a sparse logical vector.
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
#' @param direct logical; only consider if a state is reachable from 
#'   a different state (not necessarily from a start state). This is 
#'  much faster and useful for gridworlds.
#' @returns  `unreachable_states()` returns a
#'  logical vector indicating the unreachable states.
#' @export
unreachable_states <- function(model, direct = FALSE, sparse = FALSE, ...) {
  UseMethod("unreachable_states")
}


# FIXME: * Matrix::which will no longer be necessary when
#          Matrix implements the subsetting
#        * Add names to sparse vectors once matrix adds it.
#' @export
unreachable_states.MDP <- function(model,
                                   direct = FALSE,
                                   sparse = FALSE,
                                   ...) {
  
  if (!is.null(model$unreachable_states)) {
    return(.sparsify_vector(model$unreachable_states, 
                            sparse = sparse, 
                            names = model$states))
  }
  
  
  all_reachable <- function(model, sparse = FALSE) 
    .sparsify_vector(new("lsparseVector", length = length(model$states)), 
                     sparse = sparse, 
                     names = model$states)
  
  # all states are reachable from the start state
  if (all(start_vector(model) > 0))
    return(all_reachable(model, sparse))
  
  if (!direct) {
    # check by raising the transition matrix (for all actions added)
    # to the number of states
    r <- Reduce("|", transition_matrix(model, sparse = NULL))
    rr <- r
    for (i in seq_len(length(model$states) - 1L)) {
      # check if all are reachable already
      if (all(colSums(rr[Matrix::which(start_vector(model, sparse = is(rr, "sparseMatrix")) > 0), , drop = FALSE]) > 0))
        return(all_reachable(model, sparse))
      
      rr <- rr %*% r
    }
    
    rr <- colSums(rr[Matrix::which(start_vector(model, sparse = is(rr, "sparseMatrix")) > 0), , drop = FALSE]) > 0
    
    # convert to unreachable
    return (.sparsify_vector(!rr, sparse = sparse, names = model$states))
  }
  
  # direct reachability (this is a lot faster)
  r <- Reduce("|", transition_matrix(model, sparse = NULL))
  
  # self does not count
  diag(r) <- FALSE
  
  if (is(r, "sparseMatrix"))
    r <- colSums(r, sparseResult = TRUE) > 0
  else  
    r <- colSums(r) > 0
  
  # start states are also reachable
  r[start_vector(model, sparse = is(r, "sparseMatrix")) > 0] <- TRUE
  
  .sparsify_vector(!r, sparse, names = model$states)
}

#' @rdname unreachable_and_absorbing
#' @returns  `absorbing_states()` returns a logical vector indicating
#'    if the states are absorbing (terminal).
#' @export
absorbing_states <- function(model, sparse = FALSE, ...) {
  UseMethod("absorbing_states")
}

#' @export
absorbing_states.MDP <- function(model, sparse = FALSE, ...) {
  if (!is.null(model$absorbing_states)) {
    return(.sparsify_vector(model$absorbing_states, 
                           sparse = sparse, 
                           names = model$states))
  }
  
  absorbing <- rowSums(sapply(transition_matrix(model, sparse = NULL), diag)) == length(model$actions)
  
  .sparsify_vector(absorbing, sparse, names = model$states)
}

#' @rdname unreachable_and_absorbing
#' @returns `remove_unreachable_states()` returns a model with all
#'  unreachable states removed.
#' @export
remove_unreachable_states <- function(model, direct = FALSE) {
  reachable <- !unreachable_states(model, direct = direct)
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
  if(!is.null(model$absorbing_states)) 
    model$absorbing_states <- model$absorbing_states[reachable]
  if(!is.null(model$reachable_states)) 
    model$reachable_states <- model$reachable_states[reachable]
  
  # just check
  check_and_fix_MDP(model)
  model
}
