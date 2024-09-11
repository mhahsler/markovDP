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
#'   representation stored in the problem description which saves the time for conversion.
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
.translate_distribution <- function(prob,
                              model,
                              sparse = NULL) {
  ## producing the starting belief vector
  
  states <- model$states
  
  # NAs are an issue
  if (any(is.na(prob))) {
    warning("Some probabilities are NA!")
    return(prob)
  }
  
  if (is.matrix(prob) || is(prob, "dgCMatrix")) {
    if (ncol(prob) != length(states)) {
      stop("Number of column is not the number of states.")
    }
    colnames(prob) <- states
    return(.sparsify(prob, sparse = sparse))
  }
  
  # general checks for state names
  if (is.character(prob)) {
    if (any(is.na(match(
      prob, c(as.character(states), "-", "uniform")
    )))) {
      stop("Illegal probability format.\n",
           prob,
           "\nUnrecognized state name.")
    }
  }
  
  # dense probability vector
  # start: 0.3 0.1 0.0 0.2 0.5
  if (is.numeric(prob) &&
      length(prob) == length(states) &&
      round(sum(prob), 3) == 1) {
    if (!is.null(names(prob)) && !all(names(prob) == states)) {
      names(prob) <- states
    }
  }
  
  # sparse probability vector
  else if (is(prob, "sparseVector")) {
    # nothing to do
  }
  
  # start: uniform
  else if (is.character(prob) &&
           length(prob) == 1 &&
           prob[1] == "uniform") {
    prob <- rep(1 / length(states), times = length(states))
    names(prob) <- states
  }

  # start: 5
  # start include: 1 3
  # start: first-state
  # start include: first-state third state
  else if (is.character(prob) && prob[1] != "-" || 
           is.numeric(prob) && all(prob > 0)) {
    if (is.character(prob))
      prob <- match(prob, model$states)
    else
      prob <- as.integer(prob)
    
    if (length(prob) > length(states)) {
      stop("Illegal probability format.\n",
           prob,
           "\nToo many states specified.")
    }
    
    prob <- sparseVector(
      x = 1 / length(prob),
      i = prob,
      length = length(states)
    )
    
  # start exclude: -1 -3
  # start exclude: "-", "state_1", "state_2"
  } else if (is.character(prob) && prob[1] == "-" || 
             is.numeric(prob) && all(prob < 0)) {
    
    if (is.character(prob))
      prob <- match(prob[-1L], model$states)
    else
      prob <- -as.integer(prob)
    
    if (any(is.na(prob)) || 
        any(prob < 1) || 
        any(prob > length(states))) {
      stop("Illegal probability format.\n",
           prob,
           "\nState names need to exist or IDs need to be in [1, # of states].")
    }
   
    prob <- seq_along(model$states)[-prob] 
    prob <- sparseVector(
      x = 1 / length(prob),
      i = prob,
      length = length(states)
    )
  }
  
  else {
    stop("Illegal probability format.\n", prob)
  }
  
  return(.sparsify_vector(prob, sparse, names = states))
}


# make a matrix sparse if it has low density
.sparsify <- function(x,
                      sparse = TRUE) {
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


.sparsify_vector <- function(x, sparse = TRUE, names = NULL) {
  # NULL means as is but we make sure it is a vector
  if (is.null(sparse)) {
    if (is.vector(x))
      return(x)
    if (is.matrix(x))
      return(drop(x))
    
    return(as(x, "sparseVector"))
  }
  
  if (sparse)
    as(x, "sparseVector")
  else
    structure(as(x, "vector"), names = names)
}


#' @rdname accessors
#' @param start a start state description (see [MDP]). If `NULL` then the
#'   start vector is created using the start stored in the model.
#' @export
start_vector <- function(model, start = NULL, sparse = NULL) {
  if (is.null(start)) {
    start <- model$start
  }
  
  .translate_distribution(start, model = model, sparse = sparse)
}


value_matrix <-
  function(model,
           field,
           action = NULL,
           row = NULL,
           col = NULL,
           sparse = NULL,
           trans_keyword = TRUE) {
    ## action list of s x s matrices
    value <- model[[field]]
    
    # from functions
    if (is.function(value)) {
      return(function2value(model, field, value, action, row, col, sparse))
    }
    
    # from data.frame
    if (is.data.frame(value)) {
      return(df2value(model, value, action, row, col, sparse))
    }
    
    # from a list of matrices
    matrix2value(model, field, value, action, row, col, sparse, trans_keyword)
  }


function2value <- function(model, field, f, action, row, col, sparse = FALSE) {
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
  
  # do we have a 2 or argument function? (model argument does not count)
  two_args_f <- length(formals(f)) == 3L
  
  # Convert ids to names
  if (!is.null(row) && is.numeric(row))
    row <- model$states[row]
  if (!is.null(col) && is.numeric(col))
    col <- model$states[col]
  
  # single value
  if (!is.null(row) && !is.null(col)) {
    if (length(row) != 1L || length(col) != 1L)
      stop("Indices need to be single values!")
    
    if (two_args_f)
      return(unname(f(model, action, row)[col]))
    else
      return(f(model, action, row, col))
  }
  
  # rows/no cols
  if (!is.null(row) && is.null(col)) {
    if (two_args_f) {
      return(.sparsify_vector(f(model, action, row), 
                              sparse, 
                              names = model$states))
    } else{
      fv <- Vectorize(f, vectorize.args = c("end.state"))
      return(.sparsify_vector(fv(model, action, row, model$states), 
                              sparse, 
                              names = model$states))
    }
  }
  
  # no rows/cols
  if (is.null(row) && !is.null(col)) {
    fv <- Vectorize(f, vectorize.args = c("start.state"))
    if (two_args_f) {
      ### Note that the vectorized result is transposed! so col is the first index.
      return(.sparsify_vector(fv(model, action, model$states)[col , , drop = FALSE], 
                              sparse,
                              names = model$states))
    } else {
      return(.sparsify_vector(fv(model, action, model$states, col), 
                              sparse,
                              names = model$states))
    } 
  }
      
  # no rows/no cols
  if (two_args_f) {
    f_v <- Vectorize(f, vectorize.args = c("start.state"))
    o <- t(simplify2array(f_v(model, action, model$states), higher = FALSE))
  } else {
    f_v <- Vectorize(f, vectorize.args = c("start.state", "end.state"))
    o <- outer(
      model$states,
      model$states,
      FUN =
        function(r, c)
          f_v(model, action, r, c)
    )
  }
  
  # if the function returns not just 1 value!
  if (is.list(o)){
    trans <- outer(model$states, model$states, 
                   FUN = function(start, end) paste(start, "->", end))
    trans <- trans[lengths(o) != 1]
    stop("Something went wrong with the ", field, " function for action: ",
         sQuote(action), " with transition: ",
         paste(sQuote(trans), collapse = ", "))
  }
    
  #if (dim(o) != c(length(model$states), length(model$states)))
  #  stop(field, "returned an illegal vector!")
    
  dimnames(o) <- list(model$states, model$states)
  
  if (!is.null(sparse) && sparse)
    o <- .sparsify(o, sparse)
  o
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
      
      dimnames(m) <- list(model$states, model$states)
    }
    
    if (!is.character(m))
      m <- .sparsify(m, sparse)
    
    if (is.null(row) && is.null(col)) {
      return(m)
    }
    
    if (is.null(row)) {
      return(.sparsify_vector(t(m[, col, drop = FALSE]), 
                              sparse, 
                              names = model$states))
    }
    
    if (is.null(col)) {
      return(.sparsify_vector(m[row, , drop = FALSE], 
                              sparse, 
                              names = model$states))
    }
    
    return(m[row, col])
  }


df2value <-
  function(model, 
           df,
           action = NULL,
           row = NULL,
           col = NULL,
           sparse = NULL) {
    ## we use int indices here
    
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
    rows <- seq_along(model$states)
    cols <- seq_along(model$states)
    
    # return a matrix
    if (is.null(col) && is.null(row)) {
      df <-
        df[(is.na(df$action) | df$action == action), , drop = FALSE]
      
      value <- df[[4L]]
      rs <- as.integer(df[[2L]])
      cs <- as.integer(df[[3L]]) 
      
      m <- new("dgTMatrix",
               i = integer(0),
               j = integer(0), 
               x = numeric(0), 
               Dim = c(length(rows), length(cols)))
        
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
      
      return(m)
    }
    
    # TODO: Maybe sparse unless dense operation is faster
    
    # row vector
    if (is.null(col)) {
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
      
      # we default to sparse
      if (is.null(sparse))
        sparse <- TRUE
      
      return(.sparsify_vector(v, sparse, names = model$states))
    }
    
    if (is.null(row)) {
      if (is.numeric(col)) {
        col <- cols[col]
      }
      # row vector
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
    
    # value
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
