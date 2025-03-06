#' Helper Functions for Gridworld MDPs
#'
#' Helper functions for gridworld MDPs to convert between state names and
#' gridworld positions, and for visualizing policies.
#'
#' Gridworlds are implemented with state names `s(row,col)`, where
#' `row` and `col` are locations in the matrix representing the gridworld.
#' The default actions are `"up"`, `"right"`,  `"down"`, and  `"left"`.
#'
#'
#' ### Creating a Gridworld
#' `gw_init()` initializes a new gridworld creating a matrix
#' of states with the given dimensions. Other action names
#' can be specified, but they must have the same effects in the same order
#' as above. Blocked states (walls) and absorbing state can be defined.
#' This information can be used to build a custom gridworld MDP. Note that
#' blocked states are removed from the model description using
#' [remove_unreachable_states()].
#'
#' @name gridworld
#' @aliases gridworld gw
#' @family gridworld
#' @family MDP
#' @examples
#' # Defines states, actions and a transition model for a standard gridworld
#' gw <- gw_init(
#'   dim = c(7, 7),
#'   blocked_states = c("s(2,2)", "s(7,3)", "s(3,6)"),
#'   absorbing_states = "s(4,4)",
#'   state_labels = list("s(4,4)" = "Black Hole")
#' )
#'
#' str(gw)
#'
#' # display the state labels in the gridworld (states not represented in the
#' # model are shown as NA)
#' gw_matrix(gw)
#' gw_matrix(gw, what = "label")
#' gw_matrix(gw, what = "absorbing")
#' gw_matrix(gw, what = "unreachable") # these are actually missing from the model
#'
#' # a transition function for regular moves in the gridworld is provided
#' gw_transition_prob(gw, "right", "s(1,1)")
#' gw_transition_prob_end_state(gw, "right", "s(1,1)", "s(1,2)")
#'
#' # convert between state names and row/column indices
#' gw_s2rc("s(1,1)")
#' gw_rc2s(c(1, 1))
#'
#' # The information in gw can be used to build a custom MDP.
#'
#' # We modify the standard transition function so there is a 50% chance that
#' # you will get sucked into the black hole from the adjacent squares.
#' trans_black_hole <- function(model,
#'                              action,
#'                              start.state,
#'                              end.state) {
#'   # states around the black hole
#'   if (start.state %in% c(
#'     "s(3,3)", "s(3,4)", "s(3,5)", "s(4,3)", "s(4,5)",
#'     "s(5,3)", "s(5,4)", "s(5,5)"
#'   )) {
#'     if (end.state == "s(4,4)") {
#'       return(.5 + gw_transition_prob_end_state(model, action, start.state,
#'                                         end.state) * .5)
#'     } else {
#'       return(gw_transition_prob_end_state(model, action, start.state,
#'                                         end.state) * .5)
#'     }
#'   }
#'
#'   # use the standard gridworld movement
#'   gw_transition_prob_end_state(model, action, start.state, end.state)
#' }
#'
#' black_hole <- MDP(
#'   states = gw$states,
#'   actions = gw$actions,
#'   transition_prob = trans_black_hole,
#'   reward = rbind(R_(                      value = +1),
#'                  R_(end.state = "s(4,4)", value = -100),
#'                  R_(start.state = "s(4,4)", value = 0)
#'                  ),
#'   info = gw$info,
#'   name = "Black hole"
#' )
#'
#' black_hole
#' black_hole <- normalize_MDP(black_hole)
#'
#' gw_plot_transition_graph(black_hole)
#'
#' # solve the problem
#' sol <- solve_MDP(black_hole, error = 1)
#' gw_matrix(sol, what = "values")
#' gw_plot(sol)
#' # the optimal policy is to fly around, but avoid the black hole.
#'
#' # Build a Maze: The Dyna Maze from Chapter 8 in the RL book
#'
#' DynaMaze <- gw_maze_MDP(
#'   dim = c(6, 9),
#'   start = "s(3,1)",
#'   goal = "s(1,9)",
#'   walls = c(
#'     "s(2,3)", "s(3,3)", "s(4,3)",
#'     "s(5,6)",
#'     "s(1,8)", "s(2,8)", "s(3,8)"
#'   ),
#'   restart = TRUE,
#'   discount = 0.95,
#'   name = "Dyna Maze",
#' )
#' DynaMaze
#'
#' gw_matrix(DynaMaze)
#' gw_matrix(DynaMaze, what = "labels")
#'
#' gw_plot_transition_graph(DynaMaze)
#' # Note that the problems resets if the goal state would be reached.
#'
#' sol <- solve_MDP(DynaMaze, method = "LP:LP")
#'
#' gw_matrix(sol, what = "values")
#' gw_matrix(sol, what = "actions")
#' gw_plot(sol, states = TRUE)
#'
#' # check if we found a solution
#' gw_path(sol)
#'
#' # Read a maze from a text file
#' #   (X are walls, S is the start and G is the goal)
#'
#' # some examples are installed with the package
#' maze_dir <- system.file("mazes", package = "markovDP")
#' dir(maze_dir)
#'
#' file.show(file.path(maze_dir, "small_maze.txt"))
#'
#' maze <- gw_read_maze(file.path(maze_dir, "small_maze.txt"))
#' maze
#' gw_matrix(maze, what = "label")
#' gw_plot(maze)
#'
#' # Prioritized sweeping is especially effective for larger mazes.
#' sol <- solve_MDP(maze, method = "DP:GenPS")
#' sol
#'
#' gw_plot(sol)
#' gw_path(sol, horizon = 1000)
#'
#' # A maze can also be created directly from a character vector
#' maze <- gw_read_maze(
#'     textConnection(c("XXXXXX",
#'                      "XS  GX",
#'                      "XXXXXX")))
#' gw_plot(maze)
#'
#' # Create a small random maze
#' rand_maze <- gw_random_maze(dim = c(5, 5))
#' gw_plot(rand_maze)
#' @param dim vector of length two with the x and y extent of the gridworld.
#' @param actions vector with four action labels that move the agent up, right, down,
#'   and left.
#' @param blocked_states a vector with state labels for unreachable states.
#'     These states will be excluded.
#' @param absorbing_states a vector with state labels for absorbing states.
#' @param labels a list with labels for states.
#' @param state_labels a list with labels for states. The element names need
#'   to be state names.
#' @export
gw_init <-
  function(dim,
           actions = c("up", "right", "down", "left"),
           start = NULL,
           goal = NULL,
           absorbing_states = NULL,
           blocked_states = NULL,
           state_labels = list()) {
    # create states
    S <- as.vector(outer(
      seq_len(dim[1]),
      seq_len(dim[2]),
      FUN = function(x, y) {
        paste0("s(", x, ",", y, ")")
      }
    ))
    
    # TODO: remove all Matrix::which(...) once matrix implements the subsetting
    # translate start/goal
    if (!is.null(goal))
      goal <- S[Matrix::which(.translate_distribution(goal, S) > 0)]
    if (!is.null(start))
      start <- S[Matrix::which(.translate_distribution(start, S) > 0)]
    # translate unreachable
    if (!is.null(blocked_states) && length(blocked_states) > 0L) {
      blocked_states <- S[Matrix::which(.translate_distribution(blocked_states, S) > 0)]
      blocked_states <- setdiff(blocked_states, c(start, goal))
    } else
      blocke_states <- character(0)
    
    S <- setdiff(S, blocked_states)
    
    # Build reward data.frame to make unavailable action a reward of -Inf
    R <- R_(value = 0)
    
    # Add inf reward to cross boundaries
    #   # boundaries
    #   # actions = c("up", "right", "down", "left")
    #   R <- rbind(
    #     R,
    #     R_(
    #       action = actions[1],
    #       start.state = gw_rc2s(cbind(1, seq(dim[2]))),
    #       value = -Inf
    #     ),
    #     R_(
    #       action = actions[2],
    #       start.state = gw_rc2s(cbind(seq(dim[1]), dim[2])),
    #       value = -Inf
    #     ),
    #     R_(
    #       action = actions[3],
    #       start.state = gw_rc2s(cbind(dim[1], seq(dim[2]))),
    #       value = -Inf
    #     ),
    #     R_(
    #       action = actions[4],
    #       start.state = gw_rc2s(cbind(seq(dim[1]), 1)),
    #       value = -Inf
    #     )
    #   )
    # }
    
    # add start and goal to state labels
    if (!is.null(start))
      state_labels[[start]] <- state_labels[[start]] %||% "Start"
    if (!is.null(goal))
      state_labels[[goal]] <- state_labels[[goal]] %||% "Goal"
    
    l <- list(
      states = S,
      actions = actions,
      transition_prob = gw_transition_prob,
      reward = R,
      start = start,
      info = list(
        gridworld = TRUE,
        dim = dim,
        start = start,
        goal = goal,
        state_labels = state_labels
      )
    )
    
    l$absorbing_states <- normalize_state_label(absorbing_states, S)
    l$info$absorbing_states <- normalize_state_label(absorbing_states, S)
    
    l
  }



#' @rdname gridworld
#'
#' @details
#' ### Converting Between State Names and Coordinates
#'
#' `gw_s2rc()` and `gw_rc2s` help with converting from
#' state names to xy-coordinates and vice versa.
#'
#' @param model,x a solved gridworld MDP.
#' @param s a state label or a vector of labels.
#' @param rc a vector of length two with the row and column coordinate of a
#'   state in the gridworld matrix. A matrix with one state per row can be also
#'   supplied.
#' @export
gw_s2rc <- function(s) {
  if (length(s) > 1) {
    rcs <- t(sapply(s, gw_s2rc))
    rownames(rcs) <- s
    return(rcs)
  }
  
  rc <- as.integer(strsplit(s, "s\\(|,|\\)")[[1]][-1])
  if (length(rc) != 2 || any(is.na(rc))) {
    stop("Malformed gridworld state label ",
         sQuote(s),
         ". Needs to be 's(<row>,<col>)'.")
  }
  
  rc
}

#' @rdname gridworld
#' @export
gw_rc2s <- function(rc) {
  if (is.matrix(rc))
    paste0("s(", rc[, 1], ",", rc[, 2], ")")
  else
    paste0("s(", rc[1], ",", rc[2], ")")
}


#' @rdname gridworld
#'
#' @details
#' ### Inspecting Gridworlds
#' `gw_matrix()` returns different information
#'  (state names, values, actions, etc.) as a matrix. Note that some gridworlds
#'  have unreachable states removed. These states will be represented in the
#'  matrix as  `NA`.
#'
#' @param what What should be returned in the matrix. Options are:
#'  `"states"`, `"index"`, `"labels"`, `"values"`, `"actions"`, `"absorbing"`, and
#'  `"unreachable"`.
#' @export
gw_matrix <- function(model, epoch = 1L, what = "states") {
  if (is.null(model$info$gridworld)) {
    stop("'model' does not seem to be a gridworld!")
  }
  
  what <- match.arg(
    what,
    c(
      "states",
      "index",
      "labels",
      "values",
      "actions",
      "absorbing",
      "unreachable",
      "visit_probability"
    )
  )
  
  # for lists from gw_init()
  if (!inherits(model, "MDP") && !inherits(model, "MDPTF")) {
    class(model) <- "MDP"
  }
  
  nrows <- model$info$dim[1]
  ncols <- model$info$dim[2]
  all_states <- gw_init(dim = c(nrows, ncols))$states
  
  x <- switch(
    what,
    states = {
      l <- structure(rep(NA_character_, length(all_states)), names = all_states)
      l[S(model)] <- S(model)
      l
    },
    index = {
      l <- structure(rep(NA_integer_, length(all_states)), names = all_states)
      l[S(model)] <- seq_along(S(model))
      l
    },
    labels = {
      l <- structure(rep("", length(all_states)), names = all_states)
      
      # X is for states that are unreachable (incl not in the model states)
      l[!(all_states %in% S(model))] <- "X"
      labels <- model$info$state_labels
      l[names(labels)] <- unlist(labels)
      
      l
    },
    values = {
      l <- structure(rep(NA_real_, length(all_states)), names = all_states)
      p <- policy(model, drop = FALSE)[[epoch]]
      l[p$state] <- p$V
      l
    },
    actions = {
      l <- structure(rep(NA_character_, length(all_states)), names = all_states)
      p <- policy(model, drop = FALSE)[[epoch]]
      l[p$state] <- as.character(p$action)
      l
    },
    absorbing = {
      l <- structure(rep(NA, length(all_states)), names = all_states)
      l[S(model)] <- FALSE
      l[absorbing_states(model, sparse = "states")] <- TRUE
      l
    },
    unreachable = {
      l <- structure(rep(FALSE, length(all_states)), names = all_states)
      l[!(all_states %in% S(model))] <- TRUE
      l
    },
    visit_probability = {
      l <- structure(rep(NA, length(all_states)), names = all_states)
      l[S(model)] <- visit_probability(model)
      l
    }
  )
  
  matrix(x, nrow = nrows)
}

#' @rdname gridworld
#'
#' @details `gw_plot()` plots a gridworld.
#'
#' @param epoch epoch for unconverged finite-horizon solutions.
#' @param actions how to show actions. Options are:
#'  simple `"character"`, `"unicode"` arrows (needs to be supported by the used font),
#'  `"label"` of the action, and  `"none"` to suppress showing the action.
#' @param states logical; show state names.
#' @param index logical; show the state indices.
#' @param labels logical; show state labels.
#' @param impossible_actions logical; show the value and the action for absorbing states.
#' @param main logical; main title.
#' @param cex expansion factor for the action.
#' @param offset move the state labels out of the way (in fractions of a character width).
#' @param lines logical; draw lines to separate states.
#' @param contour logical; add value function contours.
#' @param col a colors for the utility values.
#' @param blocked_col a color used for blocked states. Use `NA` for no
#'   color.
#' @param ... further arguments are passed on to [graphics::image()].
#'
#' @importFrom graphics image text box abline
#' @importFrom grDevices hcl.colors
#' @export
gw_plot <-
  function(model,
           epoch = 1L,
           actions = "character",
           states = TRUE,
           index = FALSE,
           labels = TRUE,
           impossible_actions = FALSE,
           main = NULL,
           cex = 1,
           offset = .5,
           lines = TRUE,
           contour = FALSE,
           col = hcl.colors(100, "YlOrRd", rev = TRUE),
           blocked_col = "gray20",
           ...) {
    if (is.null(model$info$gridworld)) {
      stop("'model' does not seem to be a gridworld!")
    }
    
    x1 <- seq(1, model$info$dim[1])
    x2 <- seq(1, model$info$dim[2])
    
    actions <-
      match.arg(actions, c("character", "unicode", "label", "none"))
    solved <- is_solved_MDP(model)
    nrows <- model$info$dim[1]
    ncols <- model$info$dim[2]
    
    if (!solved && actions != "none") {
      actions <- "none"
    }
    
    if (!solved) {
      m <- matrix(0, nrow = nrows, ncol = ncols)
    } else {
      if (is.null(main)) {
        main <- paste("Policy:",
                      model$name,
                      paste0("(", model$solution$method, ")"))
      }
      
      m <-
        gw_matrix(model, epoch = epoch, what = "value")
      if (is.null(main)) {
        main <- model$name
      }
    }
    
    absorbing <- gw_matrix(model, what = "absorbing")
    unreachable <- gw_matrix(model, what = "unreachable")
    
    # hide unreachable values for unreachable states
    m[unreachable] <- NA
    
    # reorder for image plot
    m <- t(m)[, rev(seq_len(nrow(m))), drop = FALSE]
    
    m_plot <- m
    m_plot[m_plot == 0] <- NA
    
    # warns if we only have NAs
    suppressWarnings(image(
      x2, 
      x1,
      m_plot,
      main = main,
      axes = FALSE,
      col = col,
      ...
    ))
    
    if (contour) 
      suppressWarnings(contour(
        x2, 
        x1,
        m_plot, 
        add = TRUE)
      )
    
    # draw NAs (missing/unreachable states)
    if (!is.na(blocked_col)) {
      m <- is.na(m)
      if (any(m)) {
        m[!m] <- NA
        image(x2, x1, m, col = blocked_col, add = TRUE)
      }
    }
    
    # lines
    if (lines) {
      box()
      abline(h = x1 + .5)
      abline(v = x2 + .5)
    }
    
    g <- expand.grid(x1 = x1, x2 = x2)
    
  #  if (nrows < 2) {
  #    g$x1 <- 0
  #  }
  #  if (ncols < 2) {
  #    g$x2 <- 0
  #  }
    
    # actions
    if (actions != "none") {
      g$actions <-
        as.vector(gw_matrix(model, epoch = epoch, what = "actions"))
      
      if (!impossible_actions) {
        g$actions[absorbing | unreachable] <- NA
      }
      
      g$actions <- switch(
        actions,
        character = as.character(factor(
          g$actions,
          levels = model$actions,
          labels = c("^", ">", "v", "<")
        )),
        unicode = as.character(factor(
          g$actions,
          levels = model$actions,
          labels = c("\U2191", "\U2192", "\U2193", "\U2190")
        )),
        label = g$actions,
        none = NA
      )
     
      text(g$x2, rev(g$x1), g$actions, cex = cex)
    }
    
    if (states && !index) {
      g$state <-
        as.vector(gw_matrix(model, what = "states"))
    } else if (index && !states) {
      g$state <-
        as.vector(gw_matrix(model, what = "index"))
    } else {
      g$state <- paste0(as.vector(gw_matrix(model, what = "index")), ": ", as.vector(gw_matrix(model, what = "states")))
    }
    
    if (states || index) {
      text(
        g$x2,
        rev(g$x1),
        g$state,
        pos = 3,
        offset = offset,
        cex = .5 * cex
      )
    }
    
    if (labels) {
      g$labels <-
        as.vector(gw_matrix(model, what = "labels"))
      # hide X from drawing
      g$labels[g$labels == 'X'] <- ''
      text(
        g$x2,
        rev(g$x1),
        g$label,
        pos = 1,
        offset = offset,
        cex = .5 * cex
      )
    }
  }



#' @rdname gridworld
#'
#' @details `gw_plot_transition_graph()` plots the transition graph
#'   using the gridworld matrix as the layout.
#'
#' @param remove.loops logical; do not show transitions from a state back to itself.
#' @param vertex.color,vertex.shape,vertex.size,vertex.label,edge.arrow.size
#'  see `igraph::igraph.plotting` for details. Set `vertex.label = NULL` to show the
#'  state labels on the graph.
#' @param margin a single number specifying the margin of the plot. Can be used if the
#'   graph does not fit inside the plotting area.
#' @param main a main title for the plot. Defaults to the name of the problem.
#' @param continue logical; continue solving a solution.
#' @param ... further arguments are passed on to `igraph::plot.igraph()`.
#' @export
gw_plot_transition_graph <-
  function(x,
           remove.loops = TRUE,
           vertex.color = "gray",
           vertex.shape = "square",
           vertex.size = 10,
           vertex.label = NA,
           edge.arrow.size = .3,
           margin = .2,
           main = NULL,
           ...) {
    g <- transition_graph(x)
    
    layout <- t(sapply(x$states, gw_s2rc))[, 2:1] *
      cbind(rep(1, length(x$states)), -1)
    
    V(g)$color <- vertex.color
    V(g)$shape <- vertex.shape
    V(g)$size <- vertex.size
    
    if (remove.loops) {
      g <- igraph::simplify(g, remove.loops = TRUE)
    }
    
    asp <- x$info$dim[2] / x$info$dim[1]
    
    if (is.null(main)) {
      main <- paste("Transition Graph:", x$name)
    }
    
    plot(
      g,
      layout = norm_coords(layout, xmin = -asp, xmax = asp),
      rescale = FALSE,
      xlim = c(-asp * (1 + margin), asp * (1 + margin)),
      ylim = c(-(1 + margin), (1 + margin)),
      edge.arrow.size = edge.arrow.size,
      vertex.label = vertex.label,
      edge.label = NA,
      main = main,
      ...
    )
    
    invisible(g)
  }

#' @rdname gridworld
#'
#' @details `gw_animate()` applies algorithms from [solve_MDP()] iteration
#' by iteration and visualized the state utilities. This helps to understand
#' how the algorithms work.
#'
#' @param n number of iterations to animate.
#' @param method an MDP solution method for [solve_MDP()].
#' @param zlim limits for visualizing the state value.
#'
#' @returns `gw_animate()` returns the final solution invisibly.
#'
#' @export
gw_animate <- function(model,
                       method,
                       n,
                       zlim = NULL,
                       continue = FALSE,
                       ...) {
  if (is.null(model$info$gridworld)) {
    stop("'model' does not seem to be a gridworld!")
  }
  
  sol <- model
  
  for (i in seq(n)) {
    sol <- suppressWarnings(solve_MDP(
      sol,
      n = 1,
      method = method,
      continue = continue || (i != 1),
      ...
    ))
    
    if (!is.null(zlim)) {
      gw_plot(sol, sub = paste("Iteration", i), zlim = zlim)
    } else {
      gw_plot(sol, sub = paste("Iteration", i))
    }
  }
  
  invisible(sol)
}


#' @rdname gridworld
#'
#' @details
#'
#' ### Gridworld Transition Model
#'
#' The transition model is available in several forms:
#'
#' * `gw_transition_prob()` returns a dense vector for the action and start state.
#' * `gw_transition_prob_sparse()` returns a sparse vector for the action and start state.
#'      Note: creating sparse vectors is very expensive and should only be used
#'      for sparse models with a large state space.
#' * `gw_transition_prob_named()` returns only the non-zero probabilities as a named vector.
#' * `gw_transition_prob_end_state()` returns a single value for a given action, start and end state.
#'      Note: Using this function is very slow since it results in excessive function calls.
#'
#' @param action,start.state,end.state parameters for the transition function.
#'
#' @export
gw_transition_prob <- function(model, action, start.state) {
  S <- S(model)
  P <- setNames(numeric(length(S)), S)
  
  ai <- match(action, A(model))
  
  # stay in place for unknown actions
  if (is.na(ai)) {
    warning("Unknown action", action)
    P[start.state] <- 1
    return(P)
  }
  
  # absorbing states
  if (start.state %in% model$info$absorbing_states) {
    P[start.state] <- 1
    return(P)
  }
  
  # move
  rc <- gw_s2rc(start.state)
  rc <- switch(ai, rc + c(-1, 0), rc + c(0, +1), rc + c(+1, 0), rc + c(0, -1), )
  
  es <- gw_rc2s(rc)
  
  # stay in place if we would leave the gridworld or run into a wall
  if (!(es %in% S)) {
    es <- start.state
  }
  
  P[es] <- 1
  P
}

#' @rdname gridworld
#' @export
gw_transition_prob_sparse <- function(model, action, start.state) {
  a_i <- match(action, A(model))
  start_i = match(start.state, S(model))
  
  # stay in place for unknown actions
  if (is.na(a_i) || is.na(start_i)) {
    warning("Unknown action", action, "or start state", start.state)
    end_i <- start_i
  }
  
  # absorbing states
  else if (start.state %in% model$info$absorbing_states) {
    end_i <- start_i
  }
  
  else {
    # move
    rc <- gw_s2rc(start.state)
    rc <- switch(a_i, rc + c(-1, 0), rc + c(0, +1), rc + c(+1, 0), rc + c(0, -1), )
    
    es <- gw_rc2s(rc)
    end_i <- match(es, S(model))
    
    # stay in place if we would leave the gridworld
    if (is.na(end_i)) {
      end_i <- start_i
    }
  }
  
  return(sparseVector(
    x = 1,
    i = end_i,
    length = length(S(model))
  ))
}

#' @rdname gridworld
#' @export
gw_transition_prob_named <- function(model, action, start.state) {
  P <- structure(numeric(length(S(model))), names = S(model))
  
  ai <- match(action, A(model))
  
  # stay in place for unknown actions
  if (is.na(ai)) {
    warning("Unknown action", action)
    P[start.state] <- 1
    return(P)
  }
  
  # absorbing states
  if (start.state %in% model$info$absorbing_states) {
    P[start.state] <- 1
    return(P)
  }
  
  # move
  rc <- gw_s2rc(start.state)
  rc <- switch(ai, rc + c(-1, 0), rc + c(0, +1), rc + c(+1, 0), rc + c(0, -1), )
  
  es <- gw_rc2s(rc)
  
  # stay in place if we would leave the gridworld
  if (!(es %in% S(model))) {
    es <- start.state
  }
  
  P <- 1L
  names(P) <- es
  P
}


#' @rdname gridworld
#' @export
gw_transition_prob_end_state <- function(model, action, start.state, end.state) {
  ai <- match(action, A(model))
  
  # stay in place for unknown actions
  if (is.na(ai)) {
    warning("Unknown action", action)
    return(as.integer(end.state == start.state))
  }
  
  # stay in place for absorbing states
  absorbing_states <- model$absorbing_states
  if (!is.null(absorbing_states) &&
      start.state %in% absorbing_states) {
    return(as.integer(end.state == start.state))
  }
  
  # move
  rc <- gw_s2rc(start.state)
  rc <- switch(ai, rc + c(-1, 0), rc + c(0, +1), rc + c(+1, 0), rc + c(0, -1), )
  
  es <- gw_rc2s(rc)
  if (!(es %in% S(model))) {
    es <- start.state
  }
  as.integer(es == end.state)
}




#' @rdname gridworld
#' @details
#'
#' ### Mazes
#' `gw_maze_MDP()` helps to easily define maze-like gridworld MDPs.
#' By default, the goal state is absorbing, but with `restart = TRUE`, the
#' agent restarts the problem at the start state every time it reaches the goal
#' and receives the reward. Note that this implies that the goal state itself
#' becomes unreachable.
#'
#' @param start,goal labels for the start state and the goal state.
#' @param walls a vector with state labels for walls. Walls will
#'              become unreachable states.
#' @param goal_reward reward to transition to the goal state.
#' @param step_cost cost of each action that does not lead to the goal state.
#' @param restart logical; if `TRUE` then the problem automatically restarts when
#'      the agent reaches the goal state.
#' @param discount,horizon MDP discount factor, and horizon.
#' @param info A list with additional information. Has to contain the gridworld
#'      dimensions as element `dim` and can be created using `gw_init()`.
#' @param name a string to identify the MDP problem.
#' @param normalize logical; should the description be normalized for
#'      faster access using [normalize_MDP()].
#'
#' @returns `gw_maze_MDP()` returns an MDP object.
#' @export
gw_maze_MDP <- function(dim,
                        start,
                        goal,
                        walls = NULL,
                        actions = c("up", "right", "down", "left"),
                        goal_reward = 100,
                        step_cost = 1,
                        restart = FALSE,
                        discount = 1,
                        horizon = Inf,
                        info = NULL,
                        normalize = FALSE,
                        name = "Maze") {
  gw <-
    gw_init(
      dim,
      start = start,
      goal = goal,
      blocked_states = walls,
      absorbing_states = goal
    )
  
  if (!restart) {
    model <-
      MDP(
        states = gw$states,
        actions = gw$actions,
        transition_prob = gw$transition_prob,
        reward = rbind(
          R_(value = -step_cost),
          R_(end.state = gw$info$goal, value = goal_reward - step_cost),
          R_(start.state = gw$info$goal, value = 0) # staying in the goal is 0
        ),
        discount = discount,
        horizon = horizon,
        start = gw$start,
        info = gw$info,
        name = name
      )
    
  } else {
    # all actions in the goal redirect to a start state
    trans_restart <- function(model, action, start.state) {
      P <- structure(numeric(length(S(model))), names = S(model))
      
      if (start.state %in% model$info$goal) {
        P[model$info$start] <- 1 / length(model$info$start)
        return (P)
      }
      
      # regular move
      gw_transition_prob(model, action, start.state)
    }
    
    # note the goal state is now unreachable
    model <- MDP(
      states = gw$states,
      #actions = c(gw$actions, "restart"),
      actions = c(gw$actions),
      transition_prob = trans_restart,
      reward = rbind(
        R_(value = -step_cost),
        R_(end.state = gw$info$goal, value = goal_reward - step_cost),
        R_(start.state = gw$info$goal, value = 0) # restarting is free
      ),
      discount = discount,
      horizon = horizon,
      start = gw$start,
      info = gw$info,
      name = name
    )
  }
  
  model$absorbing_states <- gw$absorbing_states
  
  
  if (normalize) {
    model <- normalize_MDP(model)
  }
  
  model
}

#' @rdname gridworld
#' @export
gw_maze_MDPTF <- function(dim,
                          start,
                          goal,
                          walls = NULL,
                          actions = c("up", "right", "down", "left"),
                          goal_reward = 100,
                          step_cost = 1,
                          discount = 1,
                          horizon = Inf,
                          info = NULL,
                          normalize = FALSE,
                          name = "Maze") {
  if (!is.null(walls))
    walls <- normalize_state_features(walls, NULL)
  
  start <- normalize_state_features(start, NULL)
  goal <- normalize_state_features(goal, NULL)
  
  transition_func <- function(model, state, action) {
    if (absorbing_states(model, state))
      return(list(reward = 0, state_prime = state))
    
    action <- normalize_action_label(action, model)
    
    r <- -step_cost
    sp <- state + switch(
      action,
      up =   c(-1, 0),
      down = c(+1, 0),
      left = c(0, -1),
      right = c(0, +1)
    )
    
    # check bounds and walls
    if (any(sp < 1) ||
        any(sp > dim) ||
        !is.null(walls) &&
        any(apply(
          walls,
          MARGIN = 1,
          FUN = function(x)
            all(x == sp)
        ))) {
      sp <- state
      r <- -step_cost
    }
    
    # goal
    if (all(sp == goal))
      r <- goal_reward
    
    
    return(list(reward = r, state_prime = sp))
  }
  
  S <- 
  
  info <- list(
    gridworld = TRUE,
    dim = dim,
    start = features2state(start),
    goal = features2state(goal),
    state_labels = list()
  )
  
  info$state_labels[[features2state(start)]] <- info$state_labels[[features2state(start)]] %||% "Start"
  info$state_labels[[features2state(goal)]] <- info$state_labels[[features2state(goal)]] %||% "Goal"
  
  MDPTF(
    actions = actions,
    transition_func = transition_func,
    start = start,
    states = setdiff(as.vector(
      outer(
        seq_len(dim[1]),
        seq_len(dim[2]),
        FUN = function(x, y) {
          paste0("s(", x, ",", y, ")")
        }
      )
    ), normalize_state_label(walls, NULL)),
    absorbing_states = goal,
    discount = discount,
    horizon = horizon,
    info = info,
    name = name
  )
  
}


#' @rdname gridworld
#' @param wall_prob probability to make a tile a wall.
#' @export
gw_random_maze <- function(dim,
                           wall_prob = .2,
                           start = NULL,
                           goal = NULL,
                           normalize = FALSE) {
  n <- dim[1]
  m <- dim[2]
  if (is.na(m))
    m <- n
  
  if (is.null(start))
    start <- "s(1,1)"
  
  if (is.null(goal))
    goal <- paste0("s(", n, ",", m, ")")
  
  make_maze <- function() {
    # random walls (cannot be start or goal)
    if (wall_prob > 0)
      walls <- sort(sample(seq(n * m), size = ceiling(n * m * wall_prob)))
    else
      walls <- integer()
    
    maze <- gw_maze_MDP(
      dim = c(n, m),
      start = start,
      goal = goal,
      walls = walls,
      normalize = normalize,
      name = "Random Maze"
    )
    
    # random mazes may produce unreachable states and even unreachable goals
    remove_unreachable_states(maze)
  }
  
  maze <- make_maze()
  tries <- 10L
  
  while (!all(maze$info$goal %in% maze$states)) {
    if (tries < 1L)
      stop("Cannot create a valid maze after 10 tries. Please reduce wall_prob!")
    maze <- make_maze()
    tries <- tries <- -1L
  }
  
  maze
}


#' @rdname gridworld
#'
#' @details `gw_read_maze()` reads a maze in text format from a file
#' and converts it into a gridworld MDP.
#'
#' @param file filename for a maze text file.
#' @export
gw_read_maze <- function(file,
                         discount = 1,
                         restart = FALSE,
                         name = "Maze") {
  rl <- readLines(file)
  m <- do.call(rbind, strsplit(rl, split = ""))
  
  goal <- which(m == "G", arr.ind = TRUE)
  start <- which(m == "S", arr.ind = TRUE)
  walls <- which(m == "X", arr.ind = TRUE)
  
  start <- paste0("s(", paste0(start, collapse = ","), ")")
  goal <- paste0("s(", paste0(goal, collapse = ","), ")")
  walls <-
    apply(
      walls,
      MARGIN = 1,
      FUN = function(x) {
        paste0("s(", paste0(x, collapse = ","), ")")
      }
    )
  dim <- dim(m)
  
  maze <- gw_maze_MDP(
    dim = dim,
    start = start,
    goal = goal,
    walls = walls,
    restart = restart,
    discount = discount,
    name = name
  )
  maze
}


#' @rdname gridworld
#' @details `gw_path()` checks if a solved gridworld has a policy that
#' leads from the start to the goal. Note this function currently samples only a single path which is
#' an issue with stochastic transitions!
#'
#' @param start,goal start and goal states. If `NULL` then the states specified
#' in the model are used.
#'
#' @return `gw_path()` returns a list with the elements `"path"`,
#' `"reward"` and `"solved"`.
#' @export
gw_path <- function(model,
                    start = NULL,
                    goal = NULL,
                    horizon = NULL) {
  ### TODO: stochastic transitions!
  
  if (is.null(goal))
    goal <- model$info$goal
  if (is.null(goal))
    stop("No goal specified!")
  
  s <- sample_MDP(
    model,
    n = 1,
    horizon = horizon,
    epsilon = 0,
    trajectories = TRUE
  )
  
  path <- s$trajectories
  
  
  # Goal state may be non-absorbing (restart = TRUE)
  found <- match(goal, path$s_prime)
  solved <- !is.na(found)
  
  if (solved) {
    if (found == nrow(path)) {
      reward <- s$reward
    } else {
      path <- path[seq(found), , drop = FALSE]
      reward <- sum(path$r * model$discount^seq(found, 1))
    }
  } else {
    reward <- NA_real_
    path <- NULL
  }
  
  list(path = path,
       reward = reward,
       solved = solved)
}
