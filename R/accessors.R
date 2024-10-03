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
#' The accessor converts the column vectors automatically into matrices with
#' start states as rows and end states as columns. This conversion can be suppressed
#' by calling `reward_matrix(..., state_matrix = FALSE)`
#' Note that the reward structure cannot be efficiently stored using a standard sparse matrix
#' since there might be a fixed cost for each action
#' resulting in no entries with 0.
#'
#' ## Start state
#' `start_vector()` translates the start state description into a probability vector.
#'

#' @family MDP
#' @name accessors
#' @aliases accessors
#'
#' @param model A [MDP] object.
#' @param action name or index of an action.
#' @param start.state,end.state name or index of the state.
#' @param sparse logical; use sparse representation. `NULL` returns the
#'   representation stored in the problem description which saves the time
#'   for conversion. For vectors, additional options are `"states"` to return a vector of state names, and `"index"`
#'   to return an index vector of values that are `TRUE` or greater than 0.
#' @param simplify logical; try to simplify action lists into a vector or matrix?
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


# Translate a probability distribution like the start distribution in MDPs
.translate_distribution <- function(prob, state_labels, sparse = NULL, check = TRUE) {
  if (is.null(prob))
    stop("No probabilities provided (NULL)!")
  
  # NAs are an issue
  if (any(is.na(prob))) {
    warning("Some probabilities are NA!")
    return(prob)
  }
  
  if (is.matrix(prob) || is(prob, "dgCMatrix")) {
    if (ncol(prob) != length(state_labels)) {
      stop("Number of column is not the number of states.")
    }
    colnames(prob) <- state_labels
    return(.sparsify(prob, sparse = sparse))
  }
  
  # general checks for state names
  if (is.character(prob)) {
    if (any(is.na(match(
      prob, c(as.character(state_labels), "-", "uniform")
    )))) {
      stop(
        "Illegal probability format.\n",
        "[",
        paste(prob, collapse = ", "),
        "]",
        "\nUnrecognized state name."
      )
    }
  }
  
  # dense probability vector
  # start: 0.3 0.1 0.0 0.2 0.5
  if (is.numeric(prob) &&
      length(prob) == length(state_labels) &&
      sum1(prob)) {
    if (!is.null(names(prob)) && !all(names(prob) == state_labels)) {
      names(prob) <- state_labels
    }
  }
  
  # sparse or dense probability vector
  else if (is(prob, "sparseVector") || is.logical(prob)) {
    # nothing to do
  }
  
  # start: uniform
  else if (is.character(prob) &&
           length(prob) == 1 &&
           prob[1] == "uniform") {
    prob <- rep(1 / length(state_labels), times = length(state_labels))
    names(prob) <- state_labels
  }
  
  # start: 5
  # start include: 1 3
  # start: first-state
  # start include: first-state third state
  else if (is.character(prob) && (prob[1] != "-" || length(prob) == 0L) ||
           is.numeric(prob) && all(prob > 0)) {
    if (is.character(prob))
      prob <- match(prob, state_labels)
    
    if (length(prob) > length(state_labels)) {
      stop(
        "Illegal probability format.\n",
        "[",
        paste(prob, collapse = ", "),
        "]",
        "\nToo many states specified."
      )
    }
    
    i <- as.integer(prob)
    
    if (any(prob != i))
      stop(
        "Illegal probability format.\n",
        "[",
        paste(prob, collapse = ", "),
        "]",
        "\nProbabilities do not sum up to one."
      )
    
    prob <- sparseVector(
      x = rep.int(1 / length(i), length(i)),
      i = i,
      length = length(state_labels)
    )
    
    # start exclude: -1 -3
    # start exclude: "-", "state_1", "state_2"
  } else if (is.character(prob) && prob[1] == "-" ||
             is.numeric(prob) && all(prob < 0)) {
    if (is.character(prob))
      prob <- match(prob[-1L], state_labels)
    else
      prob <- -as.integer(prob)
    
    if (any(is.na(prob)) ||
        any(prob < 1) ||
        any(prob > length(state_labels))) {
      stop(
        "Illegal probability format.\n",
        "[",
        paste(prob, collapse = ", "),
        "]",
        "\nState names need to exist or IDs need to be in [1, # of states]."
      )
    }
    
    prob <- seq_along(state_labels)[-prob]
    prob <- sparseVector(
      x = 1 / length(prob),
      i = prob,
      length = length(state_labels)
    )
  }
  
  else {
    stop("Illegal probability format.\n", prob)
  }
  
  v <- .sparsify_vector(prob, sparse, names = state_labels)
  
  if (!is.character(v) && check) {
    if (!sum1(v))
      stop(
        "Illegal probability format.\n",
        "[",
        paste(prob, collapse = ", "),
        "]",
        "\nProbabilities do not sum up to one."
      )
    
    if (length(v) != length(state_labels))
      stop("Distribution does not have the correct number of entries!")
  }
  
  return(v)
}

# we could get a sparse/dense vector a set of state names, a set of state ids
.translate_logical <- function(x, state_labels, sparse = NULL) {
  # state names
  if (is.character(x)) {
    if (is.null(sparse)) {
      sparse <- "states"
    }
    
    m <- pmatch(sparse, c("states", "index")) 
    if(!is.na(m)) {
      if (m == 2L)
        x <- match(x, state_labels)
      
      return(x)
    } else {
      v <- structure(logical(length(state_labels)), names = state_labels)
      v[x] <- TRUE
      return (.sparsify_vector(v, sparse, state_labels))
    }
  }
  
  # dense/sparse vector
  if (is.logical(x) || is(x, "lsparseVector"))
    return(.sparsify_vector(x, sparse, state_labels))
  
  # index vector
  .sparsify_vector(sparseVector(rep.int(TRUE, times = length(x)), 
                                x, 
                                length = length(state_labels)), sparse, state_labels)
  
  # x <- .translate_distribution(x, state_labels, sparse, check = FALSE)
  # 
  # if(is(x, "sparseVector"))
  #   x <- as(x, "lsparseVector")
  # else if(is.numeric(x) && !is.integer(x))
  #   x <- as.logical(x)
  # 
  # x
}

# make a matrix sparse if it has low density
.sparsify <- function(x, sparse = TRUE) {
  # NULL means as is, we also keep special keywords
  if (is.null(sparse))
    return(x)
  
  if (is.character(x))
    return(x)
  
  if (!sparse) {
    if (is.matrix(x)) {
      return(x)
    } else {
      return(as.matrix(x))
    }
  }
  
  # sparse
  if (!inherits(x, "CsparseMatrix")) {
    x <- as(x, "CsparseMatrix")
  }
  
  x
}

# TODO: sparseVector currently does not have names
# accepts a sparse or dense vector
.sparsify_vector <- function(x, sparse = TRUE, names = NULL) {
  if (is.matrix(x))
    x <- drop(x)
  if (is(x, "sparseMatrix"))
    x <- as(x, "sparseVector")
  
  # NULL means as is but we make sure it is a vector
  if (is.null(sparse))
      return(x)
  
  if (is.logical(sparse)) {
    if (sparse)
      return(as(x, "sparseVector"))
    else
      return(structure(as(x, "vector"), names = names))
  }
  
  # state labels or indices
  if (is.character(sparse)) {
    sparse <- match.arg(sparse, c("states", "index"))
    if (sparse == "states") { 
      if (is(x, "sparseVector")) {
        if (is.null(names))
          stop("state names needed to return states.")
        return(names[Matrix::which(x > 0)])
      } else
        return(names(x)[x > 0])
    } else { ### index
      return(unname(Matrix::which(x > 0)))
    }
  }
  
  stop("Unknown setting for sparse.")
  
}


#' @rdname accessors
#' @param start a start state description (see [MDP]). If `NULL` then the
#'   start vector is created using the start stored in the model.
#' @export
start_vector <- function(model,
                         start = NULL,
                         sparse = NULL) {
  if (is.null(start)) {
    start <- model$start
  }
  
  if (is.null(start)) {
    start <- "uniform"
  }
  
  .translate_distribution(start, model$states, sparse = sparse)
}


value_matrix <-
  function(model,
           field,
           action = NULL,
           row = NULL,
           col = NULL,
           sparse = NULL,
           simplify = FALSE,
           trans_keyword = TRUE) {
    ## action list of s x s matrices
    value <- model[[field]]
    
    # from functions
    if (is.function(value)) {
      m <- function2value(model, field, value, action, row, col, sparse)
    }
    
    # from data.frame
    else if (is.data.frame(value)) {
      m <- df2value(model, value, action, row, col, sparse)
    }
    
    # from a list of matrices
    else {
      m <- matrix2value(model, field, value, action, row, col, sparse, trans_keyword)
    }
    
    
    if (simplify) {
      if (length(row) == 1L && length(col) == 1L)
        m <- unlist(m)
      
      else if (length(row) == 1L) {
        action_labels <- names(m)
        
        # FIXME Matrix: need to manually convert sparseVector to sparseMatrix
        if (is(m[[1]], "sparseVector"))
          m <- lapply(
            m,
            FUN = function(x)
              t(as(x, "sparseMatrix"))
          )
        
        m <- do.call(rbind, m)
        rownames(m) <- action_labels
      }
      
      else if (length(col) == 1) {
        action_labels <- names(m)
        
        # FIXME Matrix: need to manually convert sparseVector to sparseMatrix
        if (is(m[[1]], "sparseVector"))
          m <- lapply(
            m,
            FUN = function(x)
              as(x, "sparseMatrix")
          )
        
        m <- do.call(cbind, m)
        colnames(m) <- action_labels
      }
      
      else {
        # can't simplify!
      }
    }
    
    m
  }


# The user-supplied transition model function can have the
# * Arguments: model, action, start.state, end.state -> returns a single value
# * Arguments: model, action, start.state -> returns a vector of
#     length #states (possibly sparse) or
#     a short named vector with only the >0 probabilities (always dense)

function2value <- function(model,
                           field,
                           f,
                           action = NULL,
                           row = NULL,
                           col = NULL,
                           sparse = NULL) {
  
  if (is.null(action))
    action <- model$actions
  
  if (!is.character(action))
    action <- model$actions[action]
  
  if (length(action) > 1L)
    return(sapply(
      action,
      FUN = function(a)
        function2value(model, field, f, a, row, col, sparse),
      simplify = FALSE,
      USE.NAMES = TRUE
    ))
  
  # we have a single action
  
  # Convert ids to names
  if (!is.null(row) && is.numeric(row))
    row <- model$states[row]
  
  # do we have an end.state argument? (makes it 4 formal arguments)
  if (length(formals(f)) == 4L) {
    ### with end.state ####################################################
    
    if (!is.null(col) && is.numeric(col)) 
      col <- model$states[col]
    
    # single value
    if (!is.null(row) &&
        !is.null(col) &&
        length(row) == 1L &&
        length(col) == 1L) {
      return(as.numeric(f(model, action, row, col)))
    }
    
    # 1 row/no cols
    if (!is.null(row) && length(row) == 1L && is.null(col)) {
      fv <- Vectorize(f, vectorize.args = c("end.state"))
      o <- fv(model, action, row, model$states)
      return(.sparsify_vector(o, sparse = sparse, names = model$states))
    }
    
    # no rows/1 col 
    if (is.null(row) && !is.null(col) && length(col) == 1L) {
      fv <- Vectorize(f, vectorize.args = c("start.state"))
      o <- fv(model, action, model$states, col)
      return(.sparsify_vector(o, sparse = sparse, names = model$states))
    }
    
    # no rows/no cols and n rows/ n columns
    if (is.null(row))
      row <- model$states
    if (is.null(col))
      col <- model$states
    
    f_v <- Vectorize(f, vectorize.args = c("start.state", "end.state"))
    o <- outer(
      row,
      col,
      FUN =
        function(r, c)
          f_v(model, action, r, c)
    )
    
    # there is a problem if the function returns not just 1 value!
    if (is.list(o)) {
      trans <- outer(
        row,
        col,
        FUN = function(start, end)
          paste(start, "->", end)
      )
      trans <- trans[lengths(o) != 1]
      stop(
        "Something went wrong with the ",
        field,
        " function for action: ",
        sQuote(action),
        " with transition: ",
        paste(sQuote(trans), collapse = ", ")
      )
    }
    
    dimnames(o) <- list(row, col)
    
    # we default to sparse here or the matrices may become to big
    if (is.null(sparse))
      sparse <- TRUE
    
    o <- .sparsify(o, sparse)
  
    return(o)
    
  } else {
    ### no end.state ###################################################
    # the function may return a 
    # * dense probability vector 
    # * a sparseVector
    # * a short named vector with only values > 0
    # Convert it appropriately
    .f_wrapper <- function(f, model, action, row) {
      v <- f(model, action, row)
      
      # sparse or dense vector
      if(length(v) == length(model$states))
        return(v)
      
      # translate partial vector into a sparse vector
      cols <- match(names(v), model$states)
      sparseVector(x = unname(v), i = cols, length = length(model$states))
    }
    
    # we need col ids to subset sparse vectors
    if (!is.null(col) && !is.numeric(col))
      col <- match(col, model$states)
    
    # single value
    if (!is.null(row) &&
        !is.null(col) &&
        length(row) == 1L &&
        length(col) == 1L) {
      return(as.numeric(.f_wrapper(f, model, action, row)[col]))
      #return(as.numeric(f(model, action, row)[col]))
    }
    
    # 1 row /no cols
    if (!is.null(row) && length(row) == 1L &&  is.null(col)) {
      return(.sparsify_vector(.f_wrapper(f, model, action, row), sparse, names = model$states))
      #return(.sparsify_vector(f(model, action, row), sparse, names = model$states))
    }
    
    # no rows or specified rows / cols or no cols
    if (is.null(row))
      row <- model$states
    .f_wrapper_vec <- Vectorize(.f_wrapper, vectorize.args = c("row"), 
                                SIMPLIFY = FALSE)
    o <- .f_wrapper_vec(f, model, action, row)
    
    # TODO: Matrix::rbind currently does not work with sparse Vector!
    # otherwise we could use simplify = TRUE above
    if (is(o[[1]], "sparseVector"))
      o <- lapply(
        o,
        FUN = function(v)
          t(as(v, "CsparseMatrix"))
      )
    
    if (any(lengths(o) != length(model$states)))
      stop(
        "Function for ",
        field,
        " returns an incorrect number of values for action ",
        sQuote(action),
        " for start.state(s):\n\t",
        paste(sQuote(model$states[which(lengths(o) != length(model$states))]), collapse = ", ")
      )
    
    o <- do.call("rbind", o)
    dimnames(o) <- list(row, model$states)
    
    if(!is.null(col)) {
      o <- o[, col, drop = FALSE]
      colnames(o) <- model$states[col]
      
      if (length(col) == 1L)
        return(.sparsify_vector(o, sparse = sparse, names = row))
    }
    
    # we default to sparse here or the matrices may become to big
    if (is.null(sparse))
      sparse <- TRUE
    
    o <- .sparsify(o, sparse = sparse)
    
      
    return(o)
  }
}


### this just subsets the matrix list
matrix2value <-
  function(model,
           field,
           m,
           action = NULL,
           row = NULL,
           col = NULL,
           sparse = NULL,
           trans_keyword = TRUE) {
    ### TODO: It would be faster to not translate the keywords.
    if (!trans_keyword && !(is.null(row) && is.null(row)))
      trans_keyword <- TRUE
    
    if (is.null(action))
      action <- model$actions
    
    if (is.numeric(action))
      action <- model$actions[action]
    
    if (length(action) > 1L)
      return(sapply(
        action,
        FUN = function(a)
          matrix2value(model, field, m, a, row, col, sparse),
        simplify = FALSE,
        USE.NAMES = TRUE
      ))
    
    # we have a single action from here on
    m <- model[[field]][[action]]
    n_states <- length(model$states)
    
    # convert keywords
    if (trans_keyword && is.character(m)) {
      m <- switch(
        m,
        identity = {
          if (sparse %||% TRUE) {
            as(Matrix::Diagonal(n_states), "dgCMatrix")
          } else {
            diag(n_states)
          }
        },
        uniform = matrix(
          1 / length(n_states),
          nrow = length(n_states),
          ncol = length(n_states)
        )
      )
      
    }
    
    dimnames(m) <- list(model$states, model$states)
    
    ### whole matrix
    if (is.null(row) && is.null(col)) {
      return(.sparsify(m, sparse = sparse))
    }
    
    # one column
    if (is.null(row) && length(col) == 1L) {
      return(.sparsify_vector(t(m[, col, drop = FALSE]), sparse, names = model$states))
    }
    
    # one rows
    if (is.null(col) && length(row) == 1L) {
      return(.sparsify_vector(m[row, , drop = FALSE], sparse, names = model$states))
    }
    
    # single value
    if (length(col) == 1L && length(row) == 1L)
      return(m[row, col, drop = TRUE])
    
    # submatrix
    if (is.null(row))
      row <- seq_along(model$states)
    if (is.null(col))
      col <- seq_along(model$states)
    
    return(.sparsify(m[row, col, drop = FALSE], sparse = sparse))
  }


df2value <-
  function(model,
           df,
           action = NULL,
           row = NULL,
           col = NULL,
           sparse = NULL) {
    # we default to sparse
    if (is.null(sparse))
      sparse <- TRUE
    
    ## we use int indices for action here
    if (is.null(action))
      action <- model$actions
    if (is.numeric(action))
      action <- model$actions[action]
    
    if (length(action) > 1L)
      return(sapply(
        action,
        FUN = function(a)
          df2value(model, df, a, row, col, sparse),
        simplify = FALSE,
        USE.NAMES = TRUE
      ))
    
    # one action form here
    
    if (!(action %in% model$actions))
      stop("unkown action: ", action)
    
    rows <- seq_along(model$states)
    cols <- seq_along(model$states)
    
    
    # return a single value
    if (!is.null(row) &&
        !is.null(col) &&
        length(row) == 1L &&
        length(col) == 1L) {
      if (is.numeric(row)) {
        row <- rows[row]
      }
      if (is.numeric(col)) {
        col <- cols[col]
      }
      
      val <- df[[4L]][(is.na(df$action) | df$action == action) &
                        (is.na(df[[2L]]) |
                           df[[2L]] == row) &
                        (is.na(df[[3L]]) |
                           df$end.state == col)]
      
      if (length(val) == 0L) {
        return(0)
      }
      
      return(tail(val, 1L))
    }
    
    
    # TODO: Maybe sparse unless dense operation is faster
    
    # return a row vector
    if (is.null(col) && length(row) == 1L) {
      if (is.numeric(row)) {
        row <- rows[row]
      }
      df <- df[(is.na(df$action) | df$action == action) &
                 (is.na(df[[2L]]) |
                    df[[2L]] == row), , drop = FALSE]
      
      value <- df[[4L]]
      cs <- as.integer(df[[3L]])
      
      v <-
        structure(numeric(length(cols)), names = cols)
      
      for (i in seq_len(nrow(df))) {
        c <- cs[i]
        if (is.na(c)) {
          c <- cols
        }
        
        v[c] <- value[i]
      }
      

      return(.sparsify_vector(v, sparse, names = model$states))
    }
    
    # return a col vector
    if (is.null(row) && length(col) == 1L) {
      if (is.numeric(col)) {
        col <- cols[col]
      }
      df <- df[(is.na(df$action) | df$action == action) &
                 (is.na(df[[2L]]) |
                    df[[2L]] == col), , drop = FALSE]
      
      value <- df[[4L]]
      rs <- as.integer(df[[2L]])
      
      v <-
        structure(numeric(length(rows)), names = rows)
      
      for (i in seq_len(nrow(df))) {
        r <- rs[i]
        if (is.na(r)) {
          r <- rows
        }
        
        v[r] <- value[i]
      }
      
      # we default to sparse
      if (is.null(sparse))
        sparse <- TRUE
      
      return(.sparsify_vector(v, sparse, names = model$states))
    }
    
    # return the whole matrix or a submatrix
    df <-
      df[(is.na(df$action) | df$action == action), , drop = FALSE]
    
    value <- df[[4L]]
    rs <- as.integer(df[[2L]])
    cs <- as.integer(df[[3L]])
    
    ## FIXME: Matrix has an issue with subsetting dgTMatrix so I use a dgCMatrix
    ## TODO: Make this sparse
    m <- matrix(0, nrow = length(rows), ncol = length(cols))
    
    # m <- new(
    #   "dgTMatrix",
    #   i = integer(0),
    #   j = integer(0),
    #   x = numeric(0),
    #   Dim = c(length(rows), length(cols))
    # )
    #m <- as(m, "CsparseMatrix")
    
    for (i in seq_len(nrow(df))) {
      r <- rs[i]
      if (is.na(r)) {
        r <- rows
      }
      
      c <- cs[i]
      if (is.na(c)) {
        c <- cols
      }
      
      m[r, c] <- value[i]
    }
    
    m <- .sparsify(m, sparse)
    dimnames(m) <- list(model$states, model$states)
    
    if (!is.null(row))
      m <- m[row, , drop = FALSE]
    if (!is.null(col))
      m <- m[, col, drop = FALSE]
    
    return(m)
  }
