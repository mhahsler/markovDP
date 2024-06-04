#' Access to Parts of the Model Description
#'
#' Functions to provide uniform access to different parts of the MDP
#' problem description.
#'
#' Several parts of the MDP description can be defined in different ways. In particular,
#' the fields `transition_prob`, `reward`, and `start` can be defined using matrices, data frames,
#' keywords, or functions. See [MDP] for details.
#' The functions provided here, provide unified access to the data in these fields
#' to make writing code easier.
#'
#' ## Transition Probabilities \eqn{T(s'|s,a)}
#' `transition_matrix()` accesses the transition model. The complete model
#' is a list with one element for each action. Each element contains a states x states matrix
#' with \eqn{s} (`start.state`) as rows and \eqn{s'} (`end.state`) as columns.
#' Matrices with a density below 50% can be requested in sparse format
#' (as a [Matrix::dgCMatrix-class]).
#'
#' ## Reward \eqn{R(s,s',a)}
#' `reward_matrix()` accesses the reward model.
#' The preferred representation is a data.frame with the
#' columns `action`, `start.state`, `end.state`,
#' and `value`. This is a sparse representation.
#' The dense representation is a list of lists of matrices.
#' The list levels are \eqn{a} (`action`)  and \eqn{s} (`start.state`).
#' The matrices are column vectors with rows representing \eqn{s'} (`end.state`).
#' The reward structure cannot be efficiently stored using a standard sparse matrix
#' since there might be a fixed cost for each action
#' resulting in no entries with 0.
#'
#' ## Start state
#' `start_vector()` translates the start state description into a probability vector.
#'
#' ## Convert the Complete MDP Description into a consistent form
#' `normalize_MDP()` returns a new MDP definition where `transition_prob`,
#' `reward`, and `start` are normalized.
#'
#' Also, `states`, and `actions` are ordered as given in the problem
#' definition to make safe access using numerical indices possible. Normalized
#' MDP descriptions can be
#' used in custom code that expects consistently a certain format.
#' @family MDP
#' @name accessors
#' @aliases accessors
#'
#' @param x A [MDP] object.
#' @param action name or index of an action.
#' @param start.state,end.state name or index of the state.
#' @param sparse logical; use sparse matrices when the density is below 50% and keeps data.frame representation
#'  for the reward field. `NULL` returns the
#'   representation stored in the problem description which saves the time for conversion.
#' @param trans_start logical; expand the start to a probability vector?
#' @param trans_function logical; convert functions into matrices?
#' @param trans_keyword logical; convert distribution keywords (uniform and identity)
#'  in `transition_prob` matrices?
#' @param ... further arguments are passed on.
#' 
#' @return A list or a list of lists of matrices.
#' @author Michael Hahsler
#' @examples
#' data("Maze")
#' gridworld_matrix(Maze)
#'
#' # List of |A| transition matrices. One per action in the from start.states x end.states
#' Maze$transition_prob
#' transition_matrix(Maze)
#' transition_matrix(Maze, action = "up", sparse = TRUE)
#' transition_matrix(Maze,
#'   action = "up",
#'   start.state = "s(3,1)", end.state = "s(2,1)"
#' )
#'
#' # List of list of reward matrices. 1st level is action and second level is the
#' #  start state in the form of a column vector with elements for end states.
#' Maze$reward
#' reward_matrix(Maze)
#' reward_matrix(Maze, sparse = TRUE)
#' reward_matrix(Maze,
#'   action = "up",
#'   start.state = "s(3,1)", end.state = "s(2,1)"
#' )
#'
#' # Translate the initial start probability vector
#' Maze$start
#' start_vector(Maze)
#'
#' # Normalize the whole model using dense representation
#' Maze_norm <- normalize_MDP(Maze, sparse = FALSE)
#' Maze_norm$transition_prob
NULL

#' @rdname accessors
#' @param start a start state description (see [MDP]). If `NULL` then the
#'   start vector is created using the start stored in the model.
#' @export
start_vector <- function(x, start = NULL) {
  if (is.null(start)) {
    start <- x$start
  }
  .translate_belief(start, model = x)
}

# this is the start distribution in MDPs
# translate belief specifications into numeric belief vectors
.translate_belief <- function(belief = NULL, model) {
  ## producing the starting belief vector

  states <- as.character(model$states)

  if (is.null(belief)) {
    belief <- model$start
  }
  if (is.null(belief)) {
    belief <- "uniform"
  }

  if (any(is.na(belief))) {
    return(belief)
  }

  if (is.matrix(belief)) {
    if (ncol(belief) != length(states)) {
      stop("Number of column is not the number of states.")
    }
    colnames(belief) <- states
    return(belief)
  }

  # start: 0.3 0.1 0.0 0.2 0.5
  if (is.numeric(belief) &&
    length(belief) == length(states) &&
    round(sum(belief), 3) == 1) {
    if (!is.null(names(belief)) && !all(names(belief) == states)) {
      names(belief) <- states
    }
    return(belief)
  }

  # start: uniform
  if (is.character(belief) &&
    length(belief) == 1 &&
    belief[1] == "uniform") {
    belief <- rep(1 / length(states), times = length(states))
    names(belief) <- states
    return(belief)
  }


  # general checks for state IDs
  if (is.numeric(belief)) {
    belief <- as.integer(belief)
    if (any(abs(belief) < 1) || any(abs(belief) > length(states))) {
      stop(
        "Illegal belief format.\n",
        belief,
        "\nState IDs need to be in [1, # of states]."
      )
    }
  }

  # general checks for state names
  else if (is.character(belief)) {
    if (any(is.na(match(belief, c(
      as.character(states), "-"
    ))))) {
      stop(
        "Illegal belief format.\n",
        belief,
        "\nUnrecognized state name."
      )
    }
  } else {
    stop("Illegal belief format.")
  }

  # start: first-state
  # start: 5
  # start include: first-state third state
  # start include: 1 3
  if ((is.numeric(belief) && all(belief > 0)) ||
    (is.character(belief) && belief[1] != "-")) {
    if (length(belief) > length(states)) {
      stop(
        "Illegal belief format.\n",
        belief,
        "\nToo many states specified."
      )
    }
    belief_ <- rep(0, times = length(states))
    names(belief_) <- states
    belief_[belief] <- 1 / length(belief)
    return(belief_)
  }

  # start exclude: 1 3
  if (is.numeric(belief) && any(belief < 0)) {
    belief_ <- rep(1, times = length(states))
    if (length(belief) >= length(states)) {
      stop(
        "Illegal belief format.\n",
        belief,
        "\nToo many states specified."
      )
    }
    names(belief_) <- states
    belief_[-belief] <- 0
    belief_ <- belief_ / sum(belief_)
    return(belief_)
  }

  # start exclude: fifth-state seventh-state
  if (is.character(belief) && belief[1] == "-") {
    belief <- belief[-1]
    belief_ <- rep(1, times = length(states))
    names(belief_) <- states
    belief_[belief] <- 0
    belief_ <- belief_ / sum(belief_)
    return(belief_)
  }

  stop("Illegal belief format.\n", belief)
}


#' @rdname accessors
#' @export
normalize_MDP <- function(x,
                          sparse = TRUE,
                          trans_start = FALSE,
                          trans_function = TRUE,
                          trans_keyword = FALSE) {
  if (!inherits(x, "MDP")) {
    stop("x is not an MDP object!")
  }

  if (trans_start) {
    x$start <- start_vector(x)
  }

  if (is.function(x$transition_prob) && !trans_function) {
    # do nothing
  } else {
    x$transition_prob <-
      transition_matrix(x, sparse = sparse)
  }

  if ((is.function(x$reward) && !trans_function)) {
    # do nothing
  } else {
    x$reward <- reward_matrix(x, sparse = sparse)
  }

  x
}

# make a matrix sparse if it has low density
.sparsify <- function(x,
                      sparse = TRUE,
                      max_density = .5) {
  # NULL means as is, we also keep special keywords
  if (is.null(sparse) || is.character(x)) {
    return(x)
  }

  if (!sparse) {
    if (is.matrix(x)) {
      return(x)
    } else {
      return(as.matrix(x))
    }
  }

  # sparse
  if (inherits(x, "CsparseMatrix")) {
    return(x)
  }

  if (nnzero(x) / length(x) < max_density) {
    return(as(as(x, "generalMatrix"), "CsparseMatrix"))
  } else {
    return(as.matrix(x))
  }
}

value_matrix <-
  function(x,
           field,
           action = NULL,
           row = NULL,
           col = NULL,
           sparse = NULL,
           trans_keyword = TRUE) {
    ## action list of s x s matrices

    value <- x[[field]]

    # convert functions
    if (is.function(value)) {
      # shortcut for a single value
      if (!is.null(action) && !is.null(row) && !is.null(col)) {
        if (is.numeric(action)) action <- x$actions[action]
        if (is.numeric(row)) row <- x$states[row]
        if (field == "transition_prob") {
          cols <- x$states
        } else {
          ### obs
          cols <- x$observations
        }
        if (is.numeric(col)) col <- cols[col]
        return(value(action, row, col))
      }

      return(function2value(x, field, value, action, row, col, sparse))
    }

    # data.frame
    if (is.data.frame(value)) {
      return(df2value(value, action, row, col, sparse))
    }

    # we have a list of matrices
    # subset
    list2value(x, field, value, action, row, col, sparse, trans_keyword)
  }


df2value <-
  function(df,
           action = NULL,
           row = NULL,
           col = NULL,
           sparse = FALSE) {
    actions <- levels(df$action)
    rows <- levels(df[[2L]])
    cols <- levels(df[[3L]])

    if (is.null(action)) {
      l <- sapply(
        actions,
        FUN = function(a) {
          .sparsify(df2value(df, a), sparse = sparse)
        },
        simplify = FALSE
      )

      return(l)
    }

    if (is.null(col) && is.null(row)) {
      # matrix
      df <-
        df[(is.na(df$action) | df$action == action), , drop = FALSE]

      m <-
        matrix(
          0,
          nrow = length(rows),
          ncol = length(cols),
          dimnames = list(rows, cols)
        )

      for (i in seq_len(nrow(df))) {
        r <- df[[2L]][i]
        if (is.na(r)) {
          r <- rows
        }

        c <- df[[3L]][i]
        if (is.na(c)) {
          c <- cols
        }

        m[r, c] <- df$probability[i]
      }

      m <- .sparsify(m, sparse)
      return(m)
    }

    if (is.null(col)) {
      # row vector
      if (is.numeric(row)) {
        row <- rows[row]
      }
      df <- df[(is.na(df$action) | df$action == action) &
        (is.na(df[[2L]]) |
          df[[2L]] == row), , drop = FALSE]

      v <-
        structure(numeric(length(cols)), names = cols)

      for (i in seq_len(nrow(df))) {
        c <- df[[3L]][i]
        if (is.na(c)) {
          c <- cols
        }

        v[c] <- df$probability[i]
      }

      return(v)
    }

    if (is.null(row)) {
      if (is.numeric(col)) {
        col <- cols[col]
      }
      # row vector
      df <- df[(is.na(df$action) | df$action == action) &
        (is.na(df[[2L]]) |
          df[[2L]] == col), , drop = FALSE]

      v <-
        structure(numeric(length(rows)), names = rows)

      for (i in seq_len(nrow(df))) {
        r <- df[[2L]][i]
        if (is.na(r)) {
          r <- rows
        }

        v[r] <- df$probability[i]
      }

      return(v)
    }

    # value
    if (is.numeric(row)) {
      row <- rows[row]
    }
    if (is.numeric(col)) {
      col <- cols[col]
    }

    val <- df$probability[(is.na(df$action) | df$action == action) &
      (is.na(df[[2L]]) |
        df[[2L]] == row) &
      (is.na(df[[3L]]) |
        df$end.state == col)]

    if (length(val) == 0L) {
      return(0)
    }

    return(tail(val, 1L))
  }

function2value <- function(x,
                           field,
                           f,
                           action,
                           row,
                           col,
                           sparse = FALSE) {
  if (length(action) == 1L &&
    length(row) == 1L &&
    length(col) == 1L) {
    return(f(action, row, col))
  }

  # TODO: we could make access faster

  f <- Vectorize(f)
  actions <- x$actions
  rows <- x$states
  if (field == "transition_prob") {
    cols <- x$states
  } else {
    ### obs
    cols <- x$observations
  }

  m <- sapply(
    actions,
    FUN = function(a) {
      p <- outer(
        rows,
        cols,
        FUN = function(r, c) {
          f(
            a,
            r,
            c
          )
        }
      )
      dimnames(p) <- list(rows, cols)
      .sparsify(p, sparse)
    },
    simplify = FALSE
  )

  list2value(x, field, m,
    action,
    row,
    col,
    sparse = NULL
  )
}

### this just subsets the matrix list
list2value <-
  function(x,
           field,
           m,
           action = NULL,
           row = NULL,
           col = NULL,
           sparse = NULL,
           trans_keyword = TRUE) {
    actions <- x$actions
    rows <- x$states
    if (field == "transition_prob") {
      cols <- x$states
    } else {
      ### obs
      cols <- x$observations
    }

    ## convert from character
    .fix <- function(mm, sparse, trans_keyword = TRUE) {
      if (is.character(mm)) {
        if (!trans_keyword) {
          return(mm)
        }

        mm <- switch(mm,
          identity = {
            if (is.null(sparse) || sparse) {
              Matrix::Diagonal(length(rows))
            } else {
              diag(length(rows))
            }
          },
          uniform = matrix(
            1 / length(cols),
            nrow = length(rows),
            ncol = length(cols)
          )
        )

        dimnames(mm) <- list(rows, cols)
      }
      .sparsify(mm, sparse)
    }

    if (is.null(action)) {
      m <- lapply(m, .fix, sparse = sparse, trans_keyword = trans_keyword)
      return(m)
    }

    m <- .fix(m[[action]], sparse, trans_keyword)

    if (is.null(row) && is.null(col)) {
      return(m)
    }

    if (is.null(row)) {
      row <- rows
    }
    if (is.null(col)) {
      col <- cols
    }

    return(m[row, col])
  }


df2value <-
  function(df,
           action = NULL,
           row = NULL,
           col = NULL,
           sparse = FALSE) {
    actions <- levels(df$action)
    rows <- levels(df[[2L]])
    cols <- levels(df[[3L]])

    if (is.null(action)) {
      l <- sapply(
        actions,
        FUN = function(a) {
          .sparsify(df2value(df, a), sparse = sparse)
        },
        simplify = FALSE
      )

      return(l)
    }

    if (is.null(col) && is.null(row)) {
      # matrix
      df <-
        df[(is.na(df$action) | df$action == action), , drop = FALSE]

      m <-
        matrix(
          0,
          nrow = length(rows),
          ncol = length(cols),
          dimnames = list(rows, cols)
        )

      for (i in seq_len(nrow(df))) {
        r <- df[[2L]][i]
        if (is.na(r)) {
          r <- rows
        }

        c <- df[[3L]][i]
        if (is.na(c)) {
          c <- cols
        }

        m[r, c] <- df$probability[i]
      }

      m <- .sparsify(m, sparse)
      return(m)
    }

    if (is.null(col)) {
      # row vector
      if (is.numeric(row)) {
        row <- rows[row]
      }
      df <- df[(is.na(df$action) | df$action == action) &
        (is.na(df[[2L]]) |
          df[[2L]] == row), , drop = FALSE]

      v <-
        structure(numeric(length(cols)), names = cols)

      for (i in seq_len(nrow(df))) {
        c <- df[[3L]][i]
        if (is.na(c)) {
          c <- cols
        }

        v[c] <- df$probability[i]
      }

      return(v)
    }

    if (is.null(row)) {
      if (is.numeric(col)) {
        col <- cols[col]
      }
      # row vector
      df <- df[(is.na(df$action) | df$action == action) &
        (is.na(df[[2L]]) |
          df[[2L]] == col), , drop = FALSE]

      v <-
        structure(numeric(length(rows)), names = rows)

      for (i in seq_len(nrow(df))) {
        r <- df[[2L]][i]
        if (is.na(r)) {
          r <- rows
        }

        v[r] <- df$probability[i]
      }

      return(v)
    }

    # value
    if (is.numeric(row)) {
      row <- rows[row]
    }
    if (is.numeric(col)) {
      col <- cols[col]
    }

    val <- df$probability[(is.na(df$action) | df$action == action) &
      (is.na(df[[2L]]) |
        df[[2L]] == row) &
      (is.na(df[[3L]]) |
        df$end.state == col)]

    if (length(val) == 0L) {
      return(0)
    }

    return(tail(val, 1L))
  }

function2value <- function(x,
                           field,
                           f,
                           action,
                           row,
                           col,
                           sparse = FALSE) {
  if (length(action) == 1L &&
    length(row) == 1L &&
    length(col) == 1L) {
    return(f(action, row, col))
  }

  # TODO: we could make access faster

  f <- Vectorize(f)
  actions <- x$actions
  rows <- x$states
  if (field == "transition_prob") {
    cols <- x$states
  } else {
    ### obs
    cols <- x$observations
  }

  m <- sapply(
    actions,
    FUN = function(a) {
      p <- outer(
        rows,
        cols,
        FUN = function(r, c) {
          f(
            a,
            r,
            c
          )
        }
      )
      dimnames(p) <- list(rows, cols)
      .sparsify(p, sparse)
    },
    simplify = FALSE
  )

  list2value(x, field, m,
    action,
    row,
    col,
    sparse = NULL
  )
}
