#' Helper Functions for Gridworld MDPs
#'
#' Helper functions for gridworld MDPs to convert between state names and
#' gridworld positions, and for visualizing policies.
#'
#' Gridworlds are implemented with state names `s(row,col)`, where
#' `row` and `col` are locations in the matrix representing the gridworld.
#' The actions are `"up"`, `"right"`,  `"down"`, and  `"left"`.
#'
#' `gridworld_init()` initializes a new gridworld creating a matrix
#' of states with the given dimensions. Other action names
#' can be specified, but they must have the same effects in the same order
#' as above. Unreachable states (walls) and absorbing state can be defined.
#' This information can be used to build a custom gridworld MDP.
#'
#' Several helper functions are provided
#' to use states, look at the state layout, and plot policies on the
#' gridworld.
#'
#' `gridworld_maze_MDP()` helps to easily define maze-like gridworld MDPs.
#' By default, the goal state is absorbing, but with `restart = TRUE`, the
#' agent restarts the problem at the start state every time it reaches the goal
#' and receives the reward. Note that this implies that the goal state itself
#' becomes unreachable.
#'
#' `gridworld_animate()` applies algorithms from [solve_MDP()] iteration
#' by iteration and visualized the state utilities. This helps to understand
#' how the algorithms work.
#'
#' @name gridworld
#' @aliases gridworld
#' @family gridworld
#' @family MDP
#' @examples
#' # Defines states, actions and a transition model for a standard gridworld
#' gw <- gridworld_init(
#'   dim = c(7, 7),
#'   unreachable_states = c("s(2,2)", "s(7,3)", "s(3,6)"),
#'   absorbing_states = "s(4,4)",
#'   labels = list("s(4,4)" = "Black Hole")
#' )
#'
#' gw$states
#' gw$actions
#' gw$info
#'
#' # display the state labels in the gridworld
#' gridworld_matrix(gw)
#' gridworld_matrix(gw, what = "label")
#' gridworld_matrix(gw, what = "reachable")
#' gridworld_matrix(gw, what = "absorbing")
#'
#' # a transition function for regular moves in the gridworld is provided
#' gw$transition_prob("right", "s(1,1)", "s(1,2)")
#' gw$transition_prob("right", "s(2,1)", "s(2,2)") ### we cannot move into an unreachable state
#' gw$transition_prob("right", "s(2,1)", "s(2,1)") ### but the agent stays in place
#'
#' # convert between state names and row/column indices
#' gridworld_s2rc("s(1,1)")
#' gridworld_rc2s(c(1, 1))
#'
#' # The information in gw can be used to build a custom MDP.
#'
#' # We modify the standard transition function so there is a 50% chance that
#' # you will get sucked into the black hole from the adjacent squares.
#' trans_black_hole <- function(action = NA, start.state = NA, end.state = NA) {
#'   # ignore the action next to the black hole
#'   if (start.state %in% c(
#'     "s(3,3)", "s(3,4)", "s(3,5)", "s(4,3)", "s(4,5)",
#'     "s(5,3)", "s(5,4)", "s(5,5)"
#'   )) {
#'     if (end.state == "s(4,4)") {
#'       return(.5)
#'     } else {
#'       return(gw$transition_prob(action, start.state, end.state) * .5)
#'     }
#'   }
#'
#'   # use the standard gridworld movement
#'   gw$transition_prob(action, start.state, end.state)
#' }
#'
#' black_hole <- MDP(
#'   states = gw$states,
#'   actions = gw$actions,
#'   transition_prob = trans_black_hole,
#'   reward = rbind(R_(value = +1), R_(end.state = "s(4,4)", value = -100)),
#'   info = gw$info,
#'   name = "Black hole"
#' )
#'
#' black_hole
#'
#' gridworld_plot_transition_graph(black_hole)
#'
#' # solve the problem
#' sol <- solve_MDP(black_hole)
#' gridworld_matrix(sol, what = "values")
#' gridworld_plot(sol)
#' # the optimal policy is to fly around, but avoid the black hole.
#'
#' # Build a Maze: The Dyna Maze from Chapter 8 in the RL book
#'
#' DynaMaze <- gridworld_maze_MDP(
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
#' gridworld_matrix(DynaMaze)
#' gridworld_matrix(DynaMaze, what = "labels")
#'
#' gridworld_plot_transition_graph(DynaMaze)
#' # Note that the problems resets if the goal state would be reached.
#'
#' sol <- solve_MDP(DynaMaze)
#'
#' gridworld_matrix(sol, what = "values")
#' gridworld_matrix(sol, what = "actions")
#' gridworld_plot(sol)
#' gridworld_plot(sol, states = TRUE)
#'
#' # visualize the first 3 iterations of value iteration
#' gridworld_animate(DynaMaze, method = "value", n = 3)
#'
#' # Read a maze from a text file
#' #   (X are walls, S is the start and G is the goal)
#'
#' # some examples are installed with pom
#' maze_dir <- system.file("mazes", package = "markovDP")
#' dir(maze_dir)
#'
#' file.show(file.path(maze_dir, "small_maze.txt"))
#'
#' maze <- gridworld_read_maze(file.path(maze_dir, "small_maze.txt"))
#' maze
#' gridworld_plot(maze)
#' sol <- solve_MDP(maze, method = "lp", discount = 0.999)
#' sol
#'
#' gridworld_plot(sol)
#' @param dim vector of length two with the x and y extent of the gridworld.
#' @param action_labels vector with four action labels that move the agent up, right, down,
#'   and left.
#' @param unreachable_states a vector with state labels for unreachable states.
#'     These states will be excluded.
#' @param absorbing_states a vector with state labels for absorbing states.
#' @param labels a list with labels for states.
#' @param restrict_actions logical; mark actions that move outside the gridworld 
#'    or to an unreachable state as unreachable by adding a reward of `-Inf`. 
#'    If false, then the action results in the agent staying in place.
#' @export
gridworld_init <-
  function(dim,
           action_labels = c("up", "right", "down", "left"),
           unreachable_states = NULL,
           restrict_actions = FALSE,
           absorbing_states = NULL,
           labels = NULL) {
    S <- as.vector(outer(
      seq_len(dim[1]),
      seq_len(dim[2]),
      FUN = function(x, y) {
        paste0("s(", x, ",", y, ")")
      }
    ))
    
    if (!is.null(unreachable_states)) {
      S <- setdiff(S, unreachable_states)
    }
    
    T <- function(action,
                  start.state = NULL,
                  end.state = NULL) {
      ai <- pmatch(action, action_labels)
      
      # stay in place for unknown actions
      if (is.na(ai)) {
        # warning("Unknown action", action)
        return(as.integer(end.state == start.state))
      }
      
      if (!is.null(absorbing_states) &&
          start.state %in% absorbing_states) {
        return(as.integer(end.state == start.state))
      }
      
      rc <- gridworld_s2rc(start.state)
      rc <- switch(ai,
                   rc + c(-1, 0),
                   rc + c(0,+1),
                   rc + c(+1, 0),
                   rc + c(0,-1),)
      
      es <- gridworld_rc2s(rc)
      if (!(es %in% S)) {
        es <- start.state
      }
      as.integer(es == end.state)
    }
    
    #R <- rbind(R_(value = 0))
    R <- NULL
    
    if (restrict_actions) {
      # unreachable states
      for (s in unreachable_states) {
        rc <- gridworld_s2rc(s)
        dirs <- list(c(-1, 0), c(0,+1), c(+1, 0), c(0,-1))
        acts <- action_labels[c(3, 4, 1, 2)]
        
        for (i in seq_along(dirs)) {
          neighbor <- gridworld_rc2s(rc + dirs[[i]])
          R <- rbind(R,
                     R_(
                       action = acts[i],
                       start.state = neighbor,
                       value = -Inf
                     ))
        }
      }
      
      # boundaries
      # action_labels = c("up", "right", "down", "left")
      R <- rbind(
        R,
        R_(
          action = action_labels[1],
          start.state = gridworld_rc2s(cbind(1, seq(dim[2]))),
          value = -Inf
        ),
        R_(
          action = action_labels[2],
          start.state = gridworld_rc2s(cbind(seq(dim[1]), dim[2])),
          value = -Inf
        ),
        R_(
          action = action_labels[3],
          start.state = gridworld_rc2s(cbind(dim[1], seq(dim[2]))),
          value = -Inf
        ),
        R_(
          action = action_labels[4],
          start.state = gridworld_rc2s(cbind(seq(dim[1]), 1)),
          value = -Inf
        )
      )
    }
    
    
    list(
      states = S,
      actions = action_labels,
      transition_prob = T,
      reward = R,
      info = list(gridworld_dim = dim, gridworld_labels = labels)
    )
  }

#' @rdname gridworld
#' @param start,goal labels for the start state and the goal state.
#' @param walls a vector with state labels for walls. Walls will
#'              become unreachable states.
#' @param goal_reward reward to transition to the goal state.
#' @param step_cost cost of each action that does not lead to the goal state.
#' @param restart logical; if `TRUE` then the problem automatically restarts when
#'      the agent reaches the goal state.
#' @param discount,horizon MDP discount factor, and horizon.
#' @param info A list with additional information. Has to contain the gridworld
#'      dimensions as element `gridworld_dim`.
#' @param name a string to identify the MDP problem.
#' @export
gridworld_maze_MDP <- function(dim,
                               start,
                               goal,
                               walls = NULL,
                               action_labels = c("up", "right", "down", "left"),
                               goal_reward = 1,
                               step_cost = 0,
                               restart = FALSE,
                               discount = 0.9,
                               horizon = Inf,
                               info = NULL,
                               name = NA) {
  state_labels <-
    structure(list(rep_len("Goal", length(goal))), names = goal)
  
  if (start != "uniform") {
    state_labels <- c(state_labels, structure(list(rep_len(
      "Start", length(start)
    )), names = start))
  }
  
  gw <-
    gridworld_init(
      dim,
      unreachable_states = walls,
      absorbing_states = goal,
      labels = state_labels
    )
  
  if (!restart) {
    return(
      MDP(
        states = gw$states,
        actions = gw$actions,
        transition_prob = gw$transition_prob,
        reward = rbind(
          R_(value = -step_cost),
          R_(end.state = goal, value = goal_reward)
        ),
        discount = discount,
        horizon = horizon,
        start = start,
        info = gw$info,
        name = name
      )
    )
  }
  
  # with restarts. We add the action restart which is only available in the goal state.
  if (start == "uniform") {
    start <- gw$states
  }
  
  # find all states that can lead to the goal
  # v <- Vectorize(gw$transition_prob, vectorize.args = "start.state")
  # reach <- sapply(gw$actions, v, gw$states, goal)
  # reach <- reach[!rownames(reach) %in% goal, ]
  # lead.states <- names(which(rowSums(reach) > 0))
  
  # redirect to the start state
  trans_restart <- function(action = NA,
                            start.state = NA,
                            end.state = NA) {
    if (start.state %in% goal) {
      if (end.state %in% start) {
        return(1 / length(start))
      } else {
        return(0)
      }
    }
    
    # regular move
    gw$transition_prob(action, start.state, end.state)
  }
  
  
  # note the goal state is now unreachable
  MDP(
    states = gw$states,
    actions = c(gw$actions, "restart"),
    transition_prob = trans_restart,
    reward = rbind(
      R_(value = -step_cost),
      R_(action = "restart", value = -Inf),
      R_(end.state = goal, value = goal_reward),
      R_(start.state = goal, value = -Inf),
      R_(
        start.state = goal,
        action = "restart",
        value = 0
      )
    ),
    discount = discount,
    horizon = horizon,
    start = start,
    info = gw$info,
    name = name
  )
}

#' @rdname gridworld
#' @param model,x a solved gridworld MDP.
#' @param s a state label or a vector of labels.
#' @param rc a vector of length two with the row and column coordinate of a
#'   state in the gridworld matrix. A matrix with one state per row can be also 
#'   supplied.
#' @export
gridworld_s2rc <- function(s) {
  if (length(s) > 1)
    return(sapply(s, gridworld_s2rc))
  
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
gridworld_rc2s <- function(rc) {
  if (is.matrix(rc))
    paste0("s(", rc[, 1], ",", rc[, 2], ")")
  else
    paste0("s(", rc[1], ",", rc[2], ")")
}

#' @rdname gridworld
#' @param what What should be returned in the matrix. Options are:
#'  `"states"`, `"labels"`, `"values"`, `"actions"`, `"absorbing"`, and
#'  `"reachable"`.
#' @export
gridworld_matrix <- function(model,
                             epoch = 1L,
                             what = "states") {
  what <- match.arg(what,
                    c(
                      "states",
                      "labels",
                      "values",
                      "actions",
                      "absorbing",
                      "reachable"
                    ))
  
  if (is.null(model$info$gridworld_dim)) {
    stop("'model' does not seem to be a gridworld!")
  }
  
  if (!inherits(model, "MDP")) {
    class(model) <- "MDP"
  }
  
  nrows <- model$info$gridworld_dim[1]
  ncols <- model$info$gridworld_dim[2]
  
  all_states <- gridworld_init(dim = c(nrows, ncols))$states
  
  x <- switch(
    what,
    states = {
      l <- structure(rep(NA_character_, length(all_states)),
                     names = all_states)
      l[model$states] <- model$states
      l
    },
    labels = {
      l <- structure(rep("", length(all_states)),
                     names = all_states)
      
      # X is for states that are not reachable (incl not in the model states)
      l[!(all_states %in% model$states)] <- "X"
      l[!reachable_states(model)] <- "X"
      labels <- model$info$gridworld_labels
      l[names(labels)] <- unlist(labels)
      
      l
    },
    values = {
      l <- structure(rep(NA_real_, length(all_states)),
                     names = all_states)
      p <- policy(model, drop = FALSE)[[epoch]]
      l[p$state] <- p$U
      l
    },
    actions = {
      l <- structure(rep(NA_character_, length(all_states)),
                     names = all_states)
      p <- policy(model, drop = FALSE)[[epoch]]
      l[p$state] <- as.character(p$action)
      l
    },
    absorbing = {
      l <- structure(rep(TRUE, length(all_states)),
                     names = all_states)
      l[model$states] <- absorbing_states(model)
      l
    },
    reachable = {
      l <- structure(rep(FALSE, length(all_states)),
                     names = all_states)
      l[model$states] <- reachable_states(model)
      l
    }
  )
  
  matrix(x, nrow = nrows)
}

#' @rdname gridworld
#' @param epoch epoch for unconverged finite-horizon solutions.
#' @param actions how to show actions. Options are:
#'  simple `"character"`, `"unicode"` arrows (needs to be supported by the used font),
#'  `"label"` of the action, and  `"none"` to suppress showing the action.
#' @param states logical; show state names.
#' @param labels logical; show state labels.
#' @param impossible_actions logical; show the value and the action for absorbing or unreachable states.
#' @param main logical; main title.
#' @param cex expansion factor for the action.
#' @param offset move the state labels out of the way (in fractions of a character width).
#' @param lines logical; draw lines to separate states.
#' @param col a colors for the utility values.
#' @param unreachable_col a color used for unreachable states. Use `NA` for no
#'   color.
#' @param ... further arguments are passed on to [graphics::image()].
#' @importFrom graphics image text box abline
#' @importFrom grDevices hcl.colors
#' @export
gridworld_plot <-
  function(model,
           epoch = 1L,
           actions = "character",
           states = FALSE,
           labels = TRUE,
           impossible_actions = FALSE,
           main = NULL,
           cex = 1,
           offset = .5,
           lines = TRUE,
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           unreachable_col = "black",
           ...) {
    actions <-
      match.arg(actions, c("character", "unicode", "label", "none"))
    solved <- is_solved_MDP(model)
    nrows <- model$info$gridworld_dim[1]
    ncols <- model$info$gridworld_dim[2]
    
    if (!solved && actions != "none") {
      actions <- "none"
    }
    
    if (solved) {
      if (is.null(main)) {
        main <- paste("Policy:",
                      model$name,
                      paste0("(", model$solution$method, ")"))
      }
      
      U <-
        gridworld_matrix(model, epoch = epoch, what = "value")
    } else {
      if (is.null(main)) {
        main <- model$name
      }
      
      U <- matrix(NA_real_, nrow = nrows, ncol = ncols)
    }
    
    # warns if we only have NAs
    suppressWarnings(image(
      t(U)[, rev(seq_len(nrow(U))), drop = FALSE],
      main = main,
      axes = FALSE,
      col = col,
      ...
    ))
    
    # overplot unreachable squares
    if (!is.null(unreachable_col) && !is.na(unreachable_col)) {
      unreach <- matrix(NA_real_, nrow = nrows, ncol = ncols)
      unreach[!gridworld_matrix(model, what = "reachable")] <- 1
      image(t(unreach)[, rev(seq_len(nrow(unreach))), drop = FALSE],
            add = TRUE, col = unreachable_col)
    }
    
    # lines
    if (lines) {
      box()
      frac <- 1 / ((nrows - 1L) * 2)
      abline(h = (((1:nrows) * 2L) - 1L) * frac)
      frac <- 1 / ((ncols - 1L) * 2)
      abline(v = (((1:ncols) * 2L) - 1L) * frac)
    }
    
    g <- expand.grid(y = (nrows - 1L):0 / (nrows - 1L),
                     x = 0:(ncols - 1L) / (ncols - 1L))
    
    if (nrows < 2) {
      g$y <- 0
    }
    if (ncols < 2) {
      g$x <- 0
    }
    
    # states
    g$state <-
      as.vector(gridworld_matrix(model, what = "states"))
    g$labels <-
      as.vector(gridworld_matrix(model, what = "labels"))
   
    # hide X from drawing
    g$labels[g$labels == 'X'] <- ''
     
    # actions
    if (actions != "none") {
      g$actions <-
        as.vector(gridworld_matrix(model, epoch = epoch, what = "actions"))
      
      if (!impossible_actions) {
        g$actions[gridworld_matrix(model, what = "absorbing") |
                    !gridworld_matrix(model, what = "reachable")] <-
          NA
      }
      
      g$actions <- switch(
        actions,
        character = as.character(factor(
          g$actions,
          levels = c("up", "right", "down", "left"),
          labels = c("^", ">", "v", "<")
        )),
        unicode = as.character(factor(
          g$actions,
          levels = c("up", "right", "down", "left"),
          labels = c("\U2191", "\U2192", "\U2193", "\U2190")
        )),
        label = g$actions,
        none = NA
      )
      
      text(g$x, g$y, g$actions, cex = cex)
    }
    
    if (states) {
      text(
        g$x,
        g$y,
        g$state,
        pos = 3,
        offset = offset,
        cex = .5 * cex
      )
    }
    # text(g$x, g$y, g$state, pos = 3, cex = .8)
    
    if (labels) {
      text(
        g$x,
        g$y,
        g$label,
        pos = 1,
        offset = offset,
        cex = .5 * cex
      )
    }
  }

#' @rdname gridworld
#' @param hide_unreachable_states logical; do not show unreachable states.
#' @param remove.loops logical; do not show transitions from a state back to itself.
#' @param vertex.color,vertex.shape,vertex.size,vertex.label,edge.arrow.size
#'  see `igraph::igraph.plotting` for details. Set `vertex.label = NULL` to show the
#'  state labels on the graph.
#' @param margin a single number specifying the margin of the plot. Can be used if the
#'   graph does not fit inside the plotting area.
#' @param main a main title for the plot. Defaults to the name of the problem.
#' @param ... further arguments are passed on to `igraph::plot.igraph()`.
#' @export
gridworld_plot_transition_graph <-
  function(x,
           hide_unreachable_states = TRUE,
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
    
    layout <- t(sapply(x$states, gridworld_s2rc))[, 2:1] *
      cbind(rep(1, length(x$states)),-1)
    
    if (hide_unreachable_states) {
      reachable <- reachable_states(x)
      g <- induced_subgraph(g, V(g)[reachable])
      layout <- layout[reachable, , drop = FALSE]
    }
    
    V(g)$color <- vertex.color
    V(g)$shape <- vertex.shape
    V(g)$size <- vertex.size
    
    if (remove.loops) {
      g <- igraph::simplify(g, remove.loops = TRUE)
    }
    
    asp <- x$info$gridworld_dim[2] / x$info$gridworld_dim[1]
    
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
#' @param n number of iterations to animate.
#' @param method a MDP solution method for [solve_MDP()].
#' @param zlim limits for visualizing the state value.
#' @export
gridworld_animate <- function(x,
                              method,
                              n,
                              zlim = NULL,
                              ...) {
  U <- NULL
  
  for (i in seq(n)) {
    sol <- suppressWarnings(solve_MDP(
      x,
      N = 1,
      U = U,
      method = method,
      ...
    ))
    U <- value_function(sol)
    
    if (!is.null(zlim)) {
      gridworld_plot(sol,
                     sub = paste("Iteration", i),
                     zlim = zlim)
    } else {
      gridworld_plot(sol,
                     sub = paste("Iteration", i))
    }
  }
}

#' @rdname gridworld
#' @param file filename for a maze text file.
#' @export
gridworld_read_maze <- function(file,
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
  
  maze <- gridworld_maze_MDP(
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
