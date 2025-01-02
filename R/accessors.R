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
#' ## Transition Probabilities \eqn{p(s'|s,a)}
#'
#' `transition_matrix()` accesses the transition model. The complete model
#' is a list with one element for each action. Each element contains a states x states matrix
#' with \eqn{s} (`start.state`) as rows and \eqn{s'} (`end.state`) as columns.
#' Matrices with a low density can be requested in sparse format
#' (as a [Matrix::dgRMatrix-class]). It is recommended to load package `MatrixExtra`
#' to work with sparse matrices. 
#'
#' ## Reward \eqn{r(s,s',a)}
#'
#' `reward_matrix()` accesses the reward model.
#' The preferred representation is a data.frame with the
#' columns `action`, `start.state`, `end.state`,
#' and `value`. This is a sparse representation.
#'
#' The dense representation is a list of lists of matrices.
#' The list levels are \eqn{a} (`action`)  and \eqn{s} (`start.state`).
#' The matrices are column vectors with rows representing \eqn{s'} (`end.state`).
#'
#' To represent the rewards as a sparse matrices, **rewards that correspond to a transition
#' with probability zero are zeroed out if the transition model** is stored as a list
#' of matrices. This makes the reward matrices as sparse as the transition matrices.
#' The function `normalize_MDP()` with `sparse = TRUE` will perform this representation.
#'
#' ## Start state
#'
#' `start_vector()` translates the start state description into a probability vector.
#'
#' ## Sparse Matrices and Normalizing MDPs
#'
#' Different components can be specified in various ways. It is often
#' necessary to convert each component into a specific form (e.g., a
#' dense matrix) to save time during access.
#' Convert the Complete MDP Description into a consistent form
#' `normalize_MDP()` converts all components of the MDP description
#'  into a consistent form and
#' returns a new MDP definition where `transition_prob`,
#' `reward`, and `start` are normalized. This includes the internal
#' representation (dense, sparse, as a data.frame) and
#' also, `states`, and `actions` are ordered as given in the problem
#' definition to make safe access using numerical indices possible. Normalized
#' MDP descriptions can be
#' used in custom code that expects consistently a certain format.
#'
#' The default behavior of `sparse = NULL` uses parse matrices for large models
#' where the dense transition model would need more
#' than `options("MDP_SPARSE_LIMIT")` (the default is about 100 MB which
#' can be changed using
#' [options()]). Smaller models use faster dense
#' matrices.
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
#' gw_matrix(Maze)
#'
#' # here is the internal structure of the Maze object
#' str(Maze)
#'
#' # List of |A| transition matrices. One per action in the from start.states x end.states
#' Maze$transition_prob
#' transition_matrix(Maze)
#' transition_matrix(Maze, action = "up", sparse = FALSE)
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
#' start_vector(Maze, sparse = FALSE)
#' start_vector(Maze, sparse = "states")
#' start_vector(Maze, sparse = "index")
#'
#' # Normalize the whole model using sparse representation
#' Maze_norm <- normalize_MDP(Maze, sparse = TRUE)
#' str(Maze_norm)
#'
#' # Note to make the reward matrix sparse, all rewards
#' # for transitions with probability of 0 are zeroed out.
#' reward_matrix(Maze_norm)
NULL

#' @include accessors_reward.R
#' @include accessors_transitions.R

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
  
  .translate_distribution(start, S(model), sparse = sparse)
}



#' @rdname accessors
#' @param transition_prob logical; convert the transition probabilities into a list of matrices.
#' @param reward logical; convert the reward model into a list of matrices.
#' @param start logical; convert the start probability distribution into a vector.
#' @param sparse logical; use sparse matrix representation? `NULL` decides the representation
#'    based on the memory it would take to store the faster dense representation.
#' @param precompute_absorbing logical; should absorbing states be precalculated?
#' @param check_and_fix logical; checks the structure of the problem description.
#' @param progress logical; show a progress bar with estimated time for completion.
#' @export
normalize_MDP <- function(model,
                          transition_prob = TRUE,
                          reward = TRUE,
                          start = FALSE,
                          sparse = NULL,
                          precompute_absorbing = TRUE,
                          check_and_fix = FALSE,
                          progress = TRUE) {
  if (!inherits(model, "MDP")) {
    stop("model is not an MDP object!")
  }
  
  if (is.null(sparse))
    sparse <- length(S(model))^2 * length(A(model)) > getOption("MDP_SPARSE_LIMIT")
  if (!is.logical(sparse) || length(sparse) != 1L)
    stop("sparse needs to be a NULL or a logical scalar.")
  
  # start state vector + transitions matrix + reward matrix + check and fix
  n_states <- length(S(model))
  n_actions <- length(A(model))
  t_start <- n_states
  t_pass <- n_actions * n_states * n_states
  
  N <- as.numeric(start) * t_start +
    as.numeric(transition_prob) * t_pass +
    as.numeric(reward) * t_pass +
    as.numeric(precompute_absorbing &&
                 is.null(model$absorbing_states)) * t_pass +
    as.numeric(check_and_fix) * t_pass
  
  if (progress) {
    pb <- my_progress_bar(N, name = "normalize_MDP")
    pb$tick(0)
  }
  
  # start
  if (start) {
    model$start <- start_vector(model, sparse = sparse)
    if (progress)
      pb$tick(t_start)
  }
  
  # transition_prob
  if (transition_prob) {
    #model$transition_prob <- transition_matrix(model, sparse = sparse)
    #if (progress)
    #  pb$tick(t_pass)
    
    # w/progress
    model$transition_prob <-
      sapply(
        A(model),
        FUN = function(a) {
          tm <- transition_matrix(model, a, sparse = sparse)
          if (progress)
            pb$tick(n_states * n_states)
          tm
        },
        simplify = FALSE,
        USE.NAMES = TRUE
      )
  }
  
  # reward
  if (reward) {
    model$reward <- reward_matrix(model, sparse = sparse)
    if (progress)
      pb$tick(t_pass)
  }
  
  # make sure order is OK
  if (check_and_fix) {
    model <- check_and_fix_MDP(model)
    if (progress)
      pb$tick(t_pass)
  }
  
  # TODO: Normalize to state names if a different representation is used?
  # remember recalculated absorbing states
  if (precompute_absorbing && is.null(model$absorbing_states)) {
    model$absorbing_states <- absorbing_states(model, sparse = "states")
    if (progress)
      pb$tick(t_pass)
  }
  
  model
}

value_matrix <-
  function(model,
           field,
           action = NULL,
           row = NULL,
           col = NULL,
           sparse = NULL,
           drop = TRUE,
           simplify = FALSE,
           trans_keyword = TRUE) {
    if (is.null(action))
        action <- A(model)
    if(!is.character(action)) 
        action <- .normalize_action_label(action, model)
      
    # multiple actions
    if (length(action) > 1L) {
      m <- sapply(
        action,
        FUN = function(a)
          value_matrix(model, field, a, row, col, sparse, drop, simplify, trans_keyword),
        simplify = FALSE,
        USE.NAMES = TRUE
      )
    } else {
      # single action
      value <- model[[field]]
      
      # from functions
      if (is.function(value)) {
        action <- .normalize_action_label(action, model)
        row <- .normalize_state_label(row, model)
        col <- .normalize_state_label(col, model)
        m <- function2value(model, field, action, row, col, sparse, drop)
        return(m)
      }
      
      # from data.frame
      else if (is.data.frame(value)) {
        action <- .normalize_action_id(action, model)
        row <- .normalize_state_id(row, model)
        col <- .normalize_state_id(col, model)
        m <- df2value(model, field, action, row, col, sparse, drop)
        return(m)
      }
      
      # from a list of matrices
      else {
        action <- .normalize_action_id(action, model)
        row <- .normalize_state_id(row, model)
        col <- .normalize_state_id(col, model)
        m <- matrix2value(model, field, action, row, col, sparse, trans_keyword, drop)
        return(m)
      }
    }
    
    # multiple actions from here
    if (simplify && is.list(m)) {
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
      
      else if (length(col) == 1L) {
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
# Note: action, row and col need to be labels not ids!
function2value <- function(model,
                           field,
                           action,
                           row = NULL,
                           col = NULL,
                           sparse = NULL,
                           drop = TRUE) {
  
  if (length(action) != 1L)
    stop("Only single action allowed!")
  
  # we default to sparse here or the matrices may become to big
  if (is.null(sparse) &&
      length(S(model))^2 * length(A(model)) > getOption("MDP_SPARSE_LIMIT"))
    sparse <- TRUE
  
  f <- model[[field]]
  
  # Note: if we have a transition matrix with many 0s, then we can evaluate fewer
  # rewards resulting in a sparser matrix. We bypass this for single values!
  if (!(length(row) == 1L && length(col) == 1L) &&
      field == "reward" &&
      .action_is_matrix(model, "transition_prob", action)) {
    # V is a copy of T and we replace only the non-zero entries
    V <- model$transition_prob[[action]]
    if (is.null(row))
      row <- S(model)
    if (is.null(col))
      col <- S(model)
    
    V <- V[row, col, drop = FALSE]
    do <- which(V > 0, arr.ind = TRUE)
    for (i in seq_len(nrow(do)))
      V[do[i, 1L], do[i, 2L]] <- function2value(model, field, action, row[do[i, 1L]], col[do[i, 2L]])
    return(V)
  }
  
  # Convert to names
  if (!is.null(row))
    row <- as.character(row)
  
  ## do we have an end.state argument? (makes it 4 formal arguments)
  if (length(formals(f)) == 4L) {
    if (!is.null(col))
      col <- as.character(col)
    
    # single value
    if (length(row) == 1L &&
        length(col) == 1L) {
      return(as.numeric(f(model, action, row, col)))
    }
    
    # 1 row/no cols
    if (length(row) == 1L && is.null(col)) {
      f_v <- Vectorize(f, vectorize.args = c("end.state"))
      o <- f_v(model, action, row, S(model))
      if (drop)
        o <- .sparsify_vector(o, sparse = sparse, names = S(model))
      else {
        o <- rbind(o)
        rownames(o) <- row
        o <- .sparsify(o, sparse = sparse)
      }
      
      return(o)
    }
    
    # no rows/1 col
    if (is.null(row) && length(col) == 1L) {
      f_v <- Vectorize(f, vectorize.args = c("start.state"))
      o <- f_v(model, action, S(model), col)
      
      if (drop)
        o <- .sparsify_vector(o, sparse = sparse, names = S(model))
      else {
        o <- cbind(o)
        colnames(o) <- col
        o <- .sparsify(o, sparse = sparse)
      }
      
      return(o)
    }
    
    f_v <- Vectorize(f, vectorize.args = c("start.state", "end.state"))
    
    # no rows/no cols and n rows/ n columns
    if (is.null(row))
      row <- S(model)
    if (is.null(col))
      col <- S(model)
    
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
    
    o <- .sparsify(o, sparse)
    
    return(o)
    
  } else {
    ## no end.state
    # the function may return a
    # * dense probability vector
    # * a sparseVector
    # * a short named vector with only values > 0
    # we need col ids to subset sparse vectors
    if (!is.null(col) && !is.numeric(col))
      col <- match(col, S(model))
    
    # single value
    if (length(row) == 1L &&
        length(col) == 1L) {
      return(as.numeric(f(model, action, row)[col]))
    }
    
    # 1 row /no cols
    if (length(row) == 1L &&  is.null(col)) {
      o <- f(model, action, row)
      if (drop)
        o <- .sparsify_vector(o, sparse = sparse, names = S(model))
      else {
        o <- rbind(o)
        rownames(o) <- row
        o <- .sparsify(o, sparse = sparse)
      }
      
      return(o)
    }
   
    # 1 col: No advantage for function with 2 arguments
     
    # no rows or specified rows / cols or no cols
    if (is.null(row))
      row <- S(model)
    
    # Convert vector appropriately
    .f_wrapper <- function(f, model, action, row, sparse = NULL) {
      .sparsify_vector(f(model, action, row),
                       sparse = sparse,
                       names = S(model))
    }
    .f_wrapper_vec <- Vectorize(.f_wrapper,
                                vectorize.args = c("row"),
                                SIMPLIFY = FALSE)
    
    o <- .f_wrapper_vec(f, model, action, row, sparse)
     
    if (any(lengths(o) != length(S(model))))
      stop(
        "Function for ",
        field,
        " returns an incorrect number of values for action ",
        sQuote(action),
        " for start.state(s):\n\t",
        paste(sQuote(S(model)[which(lengths(o) != length(S(model)))]), collapse = ", ")
      )
    
    # Matrix/MatrixExtra has issues with rbind!
    if (is(o[[1]], "sparseVector")) {
      if (length(o) == 1L) {
        o <- o[[1]]
        
        if (drop) {
          if(!is.null(col))
            o <- o[col]
          o <- .sparsify_vector(o, sparse = sparse, names = S(model)[col])
        } else {
          o <- as.csr.matrix(o)[, col, drop = FALSE]
        }
        return(o)
      }
      #o <- Reduce(rbind2, o)
      # rbind for dsparseVectors -> dgRMatrix
      o <- do.call(MatrixExtra::rbind_csr, o)
      
      # js <- unname(lapply(o, FUN = function(x) x@i))
      # ps <- c(0L, cumsum(lengths(js)))
      # js <- unlist(js) - 1L
      # xs <- unlist(unname(lapply(o, FUN = function(x) x@x)))
      # o <- new("dgRMatrix", j = js, p = ps, x = xs, Dim = 
      #       c(length(o), length(S(model))))
    } else
    o <- do.call("rbind", o)
  
    dimnames(o) <- list(row, S(model))
    
    if (!is.null(col)) {
      o <- o[, col, drop = drop]
      #colnames(o) <- S(model)[col]
      if (length(col) == 1L)
        o <- .sparsify_vector(o, sparse = sparse, names = row)
    }
    
    return(o)
  }
}



### this just subsets the matrix list
matrix2value <-
  function(model,
           field,
           action,
           row = NULL,
           col = NULL,
           sparse = NULL,
           trans_keyword = TRUE,
           drop = TRUE) {
    ## TODO: It would be faster to not translate the keywords.
    if (!trans_keyword && !(is.null(row) && is.null(row)))
      trans_keyword <- TRUE
    
    if (length(action) != 1L)
      stop("Only single action allowed!")
    
    m <- model[[field]][[action]]
    n_states <- length(S(model))
    
    # convert keywords
    if (trans_keyword && is.character(m)) {
      m <- switch(
        m,
        identity = {
          if (sparse %||% TRUE) {
            as.csr.matrix(Matrix::Diagonal(n_states))
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
      
      dimnames(m) <- list(S(model), S(model))
      
    }
    
    ### whole matrix
    if (is.null(row) && is.null(col)) {
      return(.sparsify(m, sparse = sparse))
    }
    
    # single value
    if (length(row) == 1L && length(col) == 1L)
      return(as.numeric(m[row, col, drop = TRUE]))
    
    # one column
    if (is.null(row) && length(col) == 1L) {
      m <- m[, col, drop = drop]
      if (drop)
        m <- .sparsify_vector(m, sparse, names = S(model))
      else 
        m <- .sparsify(m , sparse)
      return(m)
    }
    
    # one row
    # TODO: Make this more memory efficient
    if (length(row) == 1L && is.null(col)) {
      m <- m[row, , drop = drop]
      if (drop)
        m <- .sparsify_vector(m, sparse, names = S(model))
      else 
        m <- .sparsify(m , sparse)
      return(m)
    }
    
    # submatrix (no drop)
    if (is.null(row))
      row <- seq_along(S(model))
    if (is.null(col))
      col <- seq_along(S(model))
    
    return(.sparsify(m[row, col, drop = FALSE], sparse = sparse))
  }


df2value <-
  function(model,
           field,
           action,
           row = NULL,
           col = NULL,
           sparse = NULL,
           drop = TRUE) {
    if (length(action) != 1L)
      stop("Only single action allowed!")
   
    df <- model[[field]]
    
    # we default to sparse for larger models
    if (is.null(sparse))
      sparse <- length(S(model))^2 * length(A(model)) > getOption("MDP_SPARSE_LIMIT")
    
    if (!is.null(row))
      row <- as.integer(row)
    if (!is.null(col))
      col <- as.integer(col)
    
    cols <- S(model)
    rows <- S(model)
    
    # return a single value fast
    if (length(row) == 1L &&
        length(col) == 1L) {
      val <- df[[4L]][(is.na(df$action) |
                         as.integer(df$action) == action) &
                        (is.na(df[[2L]]) |
                           as.integer(df[[2L]]) == row) &
                        (is.na(df[[3L]]) |
                           as.integer(df$end.state) == col)]
      
      if (length(val) == 0L) {
        return(0)
      }
      
      return(tail(val, 1L))
    }
    
    # this is slow!!!
    # Note: if we have a sparse transition matrix, then we can evaluate
    #       just P != 0
    # if (field == "reward" &&
    #     is.list(model$transition_prob)) {
    #   # we copy P and replace the non-0 entries to get V
    #   V <- model$transition_prob[[action]]
    #   # we can do sparse matrix faster ([ is slow!)
    #   if (is.matrix(V)) {
    #     if (is.null(row))
    #       row <- seq_along(S(model))
    #     if (is.null(col))
    #       col <- seq_along(S(model))
    #
    #     V <- V[row, col, drop = FALSE]
    #     do <- which(V > 0, arr.ind = TRUE)
    #     for (i in seq_len(nrow(do))) {
    #         V[do[i, 1L], do[i, 2L]] <- df2value(model, field, action, row[do[i, 1L]], col[do[i, 2L]])
    #     }
    #     return(.sparsify(V, sparse = sparse))
    #   } else if(inherits(V, "sparseMatrix")) {
    #     if (is.null(row))
    #       row <- seq_along(S(model))
    #     if (is.null(col))
    #       col <- seq_along(S(model))
    #
    #     V <- V[row, col, drop = FALSE]
    #     V <- drop0(V)
    #     #do <- which(V > 0, arr.ind = TRUE)
    #     tspm <- as(V, "TsparseMatrix")
    #     do <- cbind(tspm@i + 1L, tspm@j + 1L)
    #
    #     for (i in seq_len(nrow(do))) {
    #         #V[do[i, 1L], do[i, 2L]] <- df2value(model, field, action, row[do[i, 1L]], col[do[i, 2L]])
    #         V@x[i] <- df2value(model, field, action, row[do[i, 1L]], col[do[i, 2L]])
    #         cat("")
    #     }
    #     V <- drop0(V)
    #
    #     return(.sparsify(V, sparse = sparse))
    #   }
    # }
    
    # return a row vector
    # TODO: do this more memory efficient for sparse = TRUE!
    if (is.null(col) && length(row) == 1L) {
      df <- df[(is.na(df$action) | as.integer(df$action) == action) &
                 (is.na(df[[2L]]) |
                    as.integer(df[[2L]]) == row), , drop = FALSE]
      
      value <- df[[4L]]
      cs <- as.integer(df[[3L]])
      
      v <-
        setNames(numeric(length(cols)), cols)
      
      for (i in seq_along(cs)) {
        c <- cs[i]
        if (is.na(c)) {
          c <- cols
        }
        
        v[c] <- value[i]
      }
      
      # reward does not return values for P == 0
      if (field == "reward" &&
          .action_is_matrix(model, "transition_prob", action)) {
        v[model[["transition_prob"]][[action]][row, , drop = TRUE] == 0] <- 0
      }
      
      if (drop) {
        v <- .sparsify_vector(v, sparse, names = S(model))
      } else {
        v <- rbind(v)
        rownames(v) <- S(model)[row]
        v <- .sparsify(v, sparse)
      }
      
      return(v)
    }
    
    # return a col vector
    if (is.null(row) && length(col) == 1L) {
      df <- df[(is.na(df$action) | as.integer(df$action) == action) &
                 (is.na(df[[2L]]) |
                    as.integer(df[[2L]]) == col), , drop = FALSE]
      
      value <- df[[4L]]
      rs <- as.integer(df[[2L]])
      
      v <-
        setNames(numeric(length(rows)), rows)
      
      for (i in seq_along(rs)) {
        r <- rs[i]
        if (is.na(r)) {
          r <- rows
        }
        
        v[r] <- value[i]
      }
      
      # reward does not return values for P == 0
      if (field == "reward" &&
          .action_is_matrix(model, "transition_prob", action)) {
        v[model[["transition_prob"]][[action]][, col, drop = TRUE] == 0] <- 0
      }
      
      if (drop) {
        v <- .sparsify_vector(v, sparse, names = S(model))
      } else {
        v <- cbind(v)
        colnames(v) <- S(model)[col]
        v <- .sparsify(v, sparse)
      }
      
      return(v)
    }
    
    # return the whole matrix or a submatrix there is no more drop!
    df <-
      df[(is.na(df$action) |
            as.integer(df$action) == action), , drop = FALSE]
    
    value <- df[[4L]]
    rs <- as.integer(df[[2L]])
    cs <- as.integer(df[[3L]])
    
    ## TODO: Make this sparse for large matrices
    m <- matrix(
      0,
      nrow = length(rows),
      ncol = length(cols),
      dimnames = list(S(model), S(model))
    )
    
    # m <- new(
    #   "dgTMatrix",
    #   i = integer(0),
    #   j = integer(0),
    #   x = numeric(0),
    #   Dim = c(length(rows), length(cols))
    # )
    #m <- as(m, "CsparseMatrix")
    
    for (i in seq_along(rs)) {
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
    
    # sparse reward matrix is more efficient if we zero out entries with P of 0
    if (field == "reward" &&
        #sparse &&
        .action_is_matrix(model, "transition_prob", action)) {
      # which is for Matrix
      m[which(model[["transition_prob"]][[action]] == 0)] <- 0
    }
    
    m <- .sparsify(m, sparse)
    
    if (!is.null(row))
      m <- m[row, , drop = FALSE]
    if (!is.null(col))
      m <- m[, col, drop = FALSE]
    
    return(m)
  }
