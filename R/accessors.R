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
#' Matrices with a density below 50% can be requested in sparse format
#' (as a [Matrix::dgCMatrix-class]).
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
#' To represent the rewards as a sparse matrix, rewards that correspond to a transition
#' with probability zero are zeroed out if the transition model is stored as a list
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
    sparse <- length(S(model)) ^ 2 * length(A(model)) > getOption("MDP_SPARSE_LIMIT")
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

# Translate a probability distribution like the start distribution in MDPs
.translate_distribution <- function(prob,
                                    state_labels,
                                    sparse = NULL,
                                    check = TRUE) {
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
  else if (is.character(prob) &&
           (prob[1] != "-" || length(prob) == 0L) ||
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
  
  if (check && !is.character(v) && !is.integer(v)) {
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
    if (!is.na(m)) {
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
  .sparsify_vector(sparseVector(rep.int(TRUE, times = length(x)), x, length = length(state_labels)),
                   sparse,
                   state_labels)
  
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
  # NULL means as is
  if (is.null(sparse))
    return(x)
  
  # we also keep special keywords
  if (is.character(x))
    return(x)
  
  if (!sparse) {
    if (is.matrix(x)) {
      return(x)
    } else {
      return(as.matrix(x))
    }
  }
  
  # sparse (make it a dgCMatrix)
  if (!inherits(x, "CsparseMatrix")) {
    x <- as(as(x, "CsparseMatrix"), "generalMatrix")
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
  
  # NULL means as is
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
    } else {
      ### index
      return(unname(Matrix::which(x > 0)))
    }
  }
  
  stop("Unknown setting for sparse.")
  
}

# make sure it is a factor
.normalize_state <- function(state, model) {
  if (is.null(state) || is.factor(state))
    return(state)
  
  if (is.numeric(state)) {
    s <- factor(
      state,
      levels = seq_along(S(model)),
      labels = S(model)
    )
    if (any(is.na(s)))
      stop("Unknown state ", paste(sQuote(state[is.na(s)]), collapse = ", "))
    return(s)
  }
  
  if (is.character(state)) {
    s <- factor(state, levels = S(model))
    if (any(is.na(s)))
      stop("Unknown state ", paste(sQuote(state[is.na(s)]), collapse = ", "))
    return(s)
  } 
  stop("Unknown state label: ", sQuote(state))
}

.normalize_state_id <- function(state, model) {
  if (is.null(state) || is.integer(state))
    return(state)
  
  s <- state
  if (is.character(s)) 
    s <- match(s, S(model))
  
  if (!is.integer(s))
    s <- as.integer(s)
  
  problem <- is.na(s) | s < 1L | s > length(S(model))
  if (any(problem))
    stop("Unknown state ", paste(sQuote(state[problem]), collapse = ", "))
  
  return(s)
}

.normalize_state_label <- function(state, model) {
  if (is.null(state) || is.character(state))
    return(state)
 
  if (is.factor(state))
    return(as.character(state))
   
  return(S(model)[state])
}

.normalize_action <- function(action, model) {
  if (is.null(action) || is.factor(action))
    return(action)
 
  if (is.logical(action)) {
    if (length(action) != length(A(model)))
      stop("Illegal action definition (logical)")
    action <- which(action)
  }
   
  if (is.numeric(action)) {
    a <- factor(
      action,
      levels = seq_along(A(model)),
      labels = A(model))
    if (any(is.na(a)))
      stop("Unknown action ", paste(sQuote(action[is.na(a)]), collapse = ", "))
    return(a)
  }
  
  if (is.character(action)) {
    a <- factor(action, levels = A(model))
    if (any(is.na(a)))
      stop("Unknown action ", paste(sQuote(action[is.na(a)]), collapse = ", "))
    return(a)
  }
    
  if (is.null(action))
    return(NULL)
  
  stop("Unknown action label: ", sQuote(action))
}

.normalize_action_id <- function(action, model) {
  if (is.null(action) || is.integer(action))
    return(action)
  
  a <- action
  if (is.character(a)) 
    a <- match(a, A(model))
  
  if (!is.integer(a))
    a <- as.integer(a)
  
  problem <- is.na(a) | a < 1L | a > length(A(model))
  if (any(problem))
    stop("Unknown action ", paste(sQuote(action[problem]), collapse = ", "))
  
  return(a)
}

.normalize_action_label <- function(action, model) {
  if (is.null(action) || is.character(action))
    return(action)
  
  if (is.factor(action))
    return(as.character(action))
  
  return(A(model)[action])
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
  
  .translate_distribution(start, S(model), sparse = sparse)
}

### check if a field/action has a matrix 
.action_is_matrix <- function(model, field, action) {
  f <- model[[field]]
  
  if(is.null(f))
    stop("field ", field, " does not exist in the model!")
  
  if(is.function(f))
    return(FALSE)
  
  m <- model[[field]][[action]]
  if(is.matrix(m) || inherits(m, "Matrix"))
    return(TRUE)
  
  return(FALSE)
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
    
    value <- model[[field]]
    
    # from functions
    if (is.function(value)) {
      action <- .normalize_action_label(action, model)
      row <- .normalize_state_label(row, model)
      col <- .normalize_state_label(col, model)
      m <- function2value(model, field, action, row, col, sparse)
    }
    
    # from data.frame
    else if (is.data.frame(value)) {
      action <- .normalize_action_id(action, model)
      row <- .normalize_state_id(row, model)
      col <- .normalize_state_id(col, model)
      m <- df2value(model, field, action, row, col, sparse)
    }
    
    # from a list of matrices
    else {
      action <- .normalize_action_id(action, model)
      row <- .normalize_state_id(row, model)
      col <- .normalize_state_id(col, model)
      m <- matrix2value(model, field, action, row, col, sparse, trans_keyword)
    }
    
    
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

## TODO: For reward, we only need to evaluate the function where P !=0 !
function2value <- function(model,
                           field,
                           action = NULL,
                           row = NULL,
                           col = NULL,
                           sparse = NULL) {
  if (is.null(action))
    action <- A(model)
  
  if (length(action) > 1L)
    return(sapply(
      action,
      FUN = function(a)
        function2value(model, field, a, row, col, sparse),
      simplify = FALSE,
      USE.NAMES = TRUE
    ))
  
  # we have a single action
  f <- model[[field]]
  
  # Convert to names
  if (!is.null(row))
    row <- as.character(row)
  
  # do we have an end.state argument? (makes it 4 formal arguments)
  if (length(formals(f)) == 4L) {
    ### with end.state ####################################################
    
    if (!is.null(col))
      col <- as.character(col)
    
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
      o <- fv(model, action, row, S(model))
      return(.sparsify_vector(o, sparse = sparse, names = S(model)))
    }
    
    # no rows/1 col
    if (is.null(row) && !is.null(col) && length(col) == 1L) {
      fv <- Vectorize(f, vectorize.args = c("start.state"))
      o <- fv(model, action, S(model), col)
      return(.sparsify_vector(o, sparse = sparse, names = S(model)))
    }
    
    f_v <- Vectorize(f, vectorize.args = c("start.state", "end.state"))
    
    # no rows/no cols and n rows/ n columns
    
    # shortcut if we only want to evaluate where P != 0
    if (field == "reward" && is.null(row) && is.null(col) && .action_is_matrix(model, "transition_prob", action)) {
      m <- matrix(0, nrow = length(row), ncol = length(col), dimnames = list(row, col))
      rc <- which(model$transition_prob[[action]] != 0, arr.ind = TRUE)
      m[rc] <- f_v(model, action, rc[, 1], rc[, 2])
      
      # we default to sparse here or the matrices may become to big
      if (is.null(sparse) &&
          length(S(model)) ^ 2 * length(A(model)) > getOption("MDP_SPARSE_LIMIT"))
        sparse <- TRUE
      
      m <- .sparsify(m, sparse)
      return(m)
    }
    
    
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
    
    # we default to sparse here or the matrices may become to big
    if (is.null(sparse) &&
        length(S(model)) ^ 2 * length(A(model)) > getOption("MDP_SPARSE_LIMIT"))
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
      if (length(v) == length(S(model)))
        return(v)
      
      # translate partial vector into a sparse vector
      cols <- match(names(v), S(model))
      sparseVector(x = unname(v),
                   i = cols,
                   length = length(S(model)))
    }
    
    # we need col ids to subset sparse vectors
    if (!is.null(col) && !is.numeric(col))
      col <- match(col, S(model))
    
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
      return(.sparsify_vector(.f_wrapper(f, model, action, row), sparse, names = S(model)))
      #return(.sparsify_vector(f(model, action, row), sparse, names = S(model)))
    }
    
    # no rows or specified rows / cols or no cols
    if (is.null(row))
      row <- S(model)
    .f_wrapper_vec <- Vectorize(.f_wrapper,
                                vectorize.args = c("row"),
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
    
    if (any(lengths(o) != length(S(model))))
      stop(
        "Function for ",
        field,
        " returns an incorrect number of values for action ",
        sQuote(action),
        " for start.state(s):\n\t",
        paste(sQuote(S(model)[which(lengths(o) != length(S(model)))]), collapse = ", ")
      )
    
    o <- do.call("rbind", o)
    dimnames(o) <- list(row, S(model))
    
    if (!is.null(col)) {
      o <- o[, col, drop = FALSE]
      colnames(o) <- S(model)[col]
      
      if (length(col) == 1L)
        return(.sparsify_vector(o, sparse = sparse, names = row))
    }
    
    # we default to sparse here or the matrices may become to big
    if (is.null(sparse) &&
        length(S(model)) ^ 2 * length(A(model)) > getOption("MDP_SPARSE_LIMIT"))
      sparse <- TRUE
    
    o <- .sparsify(o, sparse = sparse)
    
    
    return(o)
  }
}



### this just subsets the matrix list
matrix2value <-
  function(model,
           field,
           action = NULL,
           row = NULL,
           col = NULL,
           sparse = NULL,
           trans_keyword = TRUE) {
    ### TODO: It would be faster to not translate the keywords.
    if (!trans_keyword && !(is.null(row) && is.null(row)))
      trans_keyword <- TRUE
    
    if (is.null(action))
      action <- seq_along(A(model))
     
    if (length(action) > 1L) {
      l <- sapply(
        action,
        FUN = function(a)
          matrix2value(model, field, a, row, col, sparse),
        simplify = FALSE
      )
      names(l) <- A(model)[action]
      return(l)
    }
    
    # we have a single action from here on
    m <- model[[field]][[action]]
    n_states <- length(S(model))
    
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
      
      dimnames(m) <- list(S(model), S(model))
      
    }
    
    ### whole matrix
    if (is.null(row) && is.null(col)) {
      return(.sparsify(m, sparse = sparse))
    }
    
    # one column
    if (is.null(row) && length(col) == 1L) {
      return(.sparsify_vector(t(m[, col, drop = FALSE]), sparse, names = S(model)))
    }
    
    # one rows
    if (is.null(col) && length(row) == 1L) {
      return(.sparsify_vector(m[row, , drop = FALSE], sparse, names = S(model)))
    }
    
    # single value
    if (length(col) == 1L && length(row) == 1L)
      return(m[row, col, drop = TRUE])
    
    # submatrix
    if (is.null(row))
      row <- seq_along(S(model))
    if (is.null(col))
      col <- seq_along(S(model))
    
    return(.sparsify(m[row, col, drop = FALSE], sparse = sparse))
  }


df2value <-
  function(model,
           field,
           action = NULL,
           row = NULL,
           col = NULL,
           sparse = NULL) {
    ## we use int indices for action here
    if (is.null(action))
      action <- seq_along(A(model))
    
    if (length(action) > 1L) {
      r <- sapply(
        action,
        FUN = function(a)
          df2value(model, field, a, row, col, sparse),
        simplify = FALSE
        )
      names(r) <- A(model)[action]
      
      return(r)
    }
    
    # one action form here
    df <- model[[field]]
     
    # we default to sparse
    if (is.null(sparse))
      sparse <- length(S(model)) ^ 2 * length(A(model)) > getOption("MDP_SPARSE_LIMIT")
    
    if (!is.null(row))
      row <- as.integer(row)
    if (!is.null(col))
      col <- as.integer(col)
    
    cols <- S(model)
    rows <- S(model)
    
    # return a single value fast
    if (!is.null(row) &&
        !is.null(col) &&
        length(row) == 1L &&
        length(col) == 1L) {
      
      val <- df[[4L]][(is.na(df$action) | as.integer(df$action) == action) &
                        (is.na(df[[2L]]) |
                           as.integer(df[[2L]]) == row) &
                        (is.na(df[[3L]]) |
                           as.integer(df$end.state) == col)]
      
      if (length(val) == 0L) {
        return(0)
      }
      
      return(tail(val, 1L))
    }
    
    
    # TODO: Maybe sparse unless dense operation is faster
    
    # return a row vector
    if (is.null(col) && length(row) == 1L) {
      df <- df[(is.na(df$action) | as.integer(df$action) == action) &
                 (is.na(df[[2L]]) |
                    as.integer(df[[2L]]) == row), , drop = FALSE]
      
      value <- df[[4L]]
      cs <- as.integer(df[[3L]])
      
      v <-
        structure(numeric(length(cols)), names = cols)
      
      for (i in seq_along(cs)) {
        c <- cs[i]
        if (is.na(c)) {
          c <- cols
        }
        
        v[c] <- value[i]
      }
      
      return(.sparsify_vector(v, sparse, names = S(model)))
    }
    
    # return a col vector
    if (is.null(row) && length(col) == 1L) {
      df <- df[(is.na(df$action) | as.integer(df$action) == action) &
                 (is.na(df[[2L]]) |
                    as.integer(df[[2L]]) == col), , drop = FALSE]
      
      value <- df[[4L]]
      rs <- as.integer(df[[2L]])
      
      v <-
        structure(numeric(length(rows)), names = rows)
      
      for (i in seq_along(rs)) {
        r <- rs[i]
        if (is.na(r)) {
          r <- rows
        }
        
        v[r] <- value[i]
      }
      
      return(.sparsify_vector(v, sparse, names = S(model)))
    }
    
    # return the whole matrix or a submatrix
    df <-
      df[(is.na(df$action) | as.integer(df$action) == action), , drop = FALSE]
    
    value <- df[[4L]]
    rs <- as.integer(df[[2L]])
    cs <- as.integer(df[[3L]])
    
    ## TODO: Make this sparse for large matrices
    m <- matrix(0, nrow = length(rows), ncol = length(cols), 
                dimnames = list(S(model), S(model)))
    
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
    if (field == "reward" && sparse && .action_is_matrix(model, "transition_prob", action)) {
      m[Matrix::which(model[["transition_prob"]][[action]] == 0)] <- 0
    }
    
    m <- .sparsify(m, sparse)
  
    if (!is.null(row))
      m <- m[row, , drop = FALSE]
    if (!is.null(col))
      m <- m[, col, drop = FALSE]
    
    return(m)
  }
