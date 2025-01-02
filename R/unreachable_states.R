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
#' @importFrom fastmap fastmap faststack
#' @export
unreachable_states <- function(model,
                               horizon = Inf,
                               sparse = "states",
                               progress = TRUE,
                               ...) {
  UseMethod("unreachable_states")
}

#' @export
unreachable_states.MDP <- function(model,
                                   horizon = Inf,
                                   sparse = "states",
                                   progress = TRUE,
                                   ...) {
  if (is.null(sparse))
    sparse <- TRUE
  
  .translate_logical(!.reachable_states(model, horizon = horizon, progress = progress),
                     S(model),
                     sparse = sparse)
}

# depth-first search
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
    pb <- my_progress_bar(N = length(S(model)), name = "unreachable_states")
    pb$tick(0)
  }
  
  while (frontier$size() > 0) {
    state <- frontier$pop()
    
    if (progress) {
      pb$tick()
    }
    
    # available_actions is slow!
    #for (action in available_actions(model, state)){
    for (action in A(model)) {
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
  
  S(model) %in% names(reached$as_list())
}

#' @rdname unreachable_states
#' @returns `remove_unreachable_states()` returns a model with all
#'  unreachable states removed.
#' @export
remove_unreachable_states <- function(model, ...) {
  # Use dense because we need a named vector!
  reachable <-  which(!unreachable_states(model, 
                                    sparse = FALSE))
  
  if (length(reachable) == length(S(model))) {
    model$unreachable <- character(0)
    return(model)
  }
  
  # remove states from a list of matrices or a data.frame
  keep_states <- function(field, states) {
    if (is.data.frame(field)) {
      field <-
        field[field$start.state %in% c(NA, names(states)) |
                field$end.state %in% c(NA, names(states)), , drop = FALSE]
      field$start.state <-
        factor(as.character(field$start.state), levels = names(states))
      field$end.state <-
        factor(as.character(field$end.state), levels = names(states))
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
    if (length(model$start) == length(S(model))) {
      ### prob vector
      model$start <- model$start[reachable]
      if (!sum1(model$start)) {
        stop(
          "Probabilities for reachable states do not sum up to one! An unreachable state had a non-zero probability."
        )
      }
    } else {
      ### state ids... we translate to state names and use code below!
      model$start <- S(model)[model$start]
    }
  }
  if (is.character(model$start)) {
    if (model$start == "uniform") {
      # do nothing
    } else {
      model$start <- intersect(model$start, S(model)[reachable])
    }
    if (length(model$start) == 0L) {
      stop("Start state is not reachable.")
    }
  }
  
  model$states <- model$states[reachable]
  model$transition_prob <- keep_states(model$transition_prob, reachable)
  model$reward <- keep_states(model$reward, reachable)
  
  # update reachable and unreachable
  if (!is.null(model$absorbing_states))
    model$absorbing_states <- interaction(absorbing_states(model, sparse = "states"), names(reachable))
  
  # just check
  check_and_fix_MDP(model)
  model
}
