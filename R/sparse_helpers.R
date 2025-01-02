
# make a matrix sparse if it has low density
.sparsify <- function(x, sparse = TRUE) {
  # NULL means as is
  if (is.null(sparse))
    return(x)
  
  # we also keep special keywords
  if (is.character(x))
    return(x)
  
  if (!sparse)
      return(as.matrix(x))

  
  x <- as.csr.matrix(x)
  x
}

# accepts a sparse or dense vector or a matrix with one row or column
.sparsify_vector <- function(x, sparse = TRUE, names = NULL) {
  if (is.matrix(x))
    x <- drop(x)
  
  if (is(x, "sparseMatrix"))
    x <- as.sparse.vector(x)
  
  # expand states format into vector
  if (!is.null(names) && length(names) != length(x))
    x <- sparseVector(
      x = unname(x),
      i = match(names(x), names),
      length = length(names)
    )
  
  # NULL means as is
  if (is.null(sparse))
    return(x)
  
  if (is.logical(sparse)) {
    if (sparse)
      return(as.sparse.vector(x))
    else
      return(setNames(as(x, "vector"), names))
  }
  
  # state labels or indices
  if (is.character(sparse)) {
    sparse <- match.arg(sparse, c("states", "index"))
    if (sparse == "states") {
      if (is(x, "sparseVector")) {
        if (is.null(names))
          stop("state names needed to return states.")
        # faster
        #return(names[which(x > 0)])
        return(names[x@i])
      } else
        return(names(x)[x > 0])
    } else {
      ### index
      return(unname(which(x > 0)))
    }
  }
  
  stop("Unknown setting for sparse.")
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
  
  if (is.matrix(prob) || is(prob, "sparseMatrix")) {
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

# translate any representation into an logical representation
# x can be:
# * a character vector with state names
# * a integer vector with state ids
# * a logical vector
# * a lsparseVector
#
# output:
# * sparse = NULL -> as is
# * sparse = TRUE -> lsparseVector
# * sparse = FALSE -> logical vector
# * sparse = "states" -> character vector 
# * sparse = "index" -> integer vector
#
.translate_logical <- function(x, state_labels, sparse = NULL) {
  if (is.matrix(x))
    x <- drop(x)
  
  if (is(x, "sparseMatrix"))
    x <- as.sparse.vector(x, logical = TRUE)
  
  if (is.character(sparse)) {
    sparse <- match.arg(sparse, c("states", "index"))
  }
  
  # NULL means as is
  if (is.null(sparse))
    return(x)
  
  if (is.character(x)) {
    if (sparse == "states")
      return(x)
    x <- sparseVector(
      x = rep.int(TRUE, length(x)),
      i = match(x, state_labels),
      length = length(state_labels)
    )
  }
  
  if (is.integer(x) || is.numeric(x)) {
    if (sparse == "index")
      return(x)
    x <- sparseVector(
      x = rep.int(TRUE, length(x)),
      i = as.integer(x),
      length = length(state_labels)
    )
  }
  
  if (is.logical(sparse)) {
    if (sparse)
      return(as.sparse.vector(x, logical = TRUE))
    else
      return(setNames(as(x, "logical"), state_labels))
  }
  
  # state labels or indices
  if (sparse == "states") {
    if (is(x, "sparseVector")) {
      if (is.null(state_labels))
        stop("state names needed to return states.")
      # faster
      #return(names[which(x > 0)])
      return(state_labels[x@i])
    } else
      return(state_labels[x > 0])
  }
  
  if (sparse == "index")
    return(unname(which(x > 0)))
  
  stop("Unable to translate vector to logical!")
}


# make sure state and actions are factors, integer vectors or character vectors
.normalize_state <- function(state, model) {
  if (is.null(state) || is.factor(state))
    return(state)
  
  if (is.numeric(state)) {
    s <- factor(state, levels = seq_along(S(model)), labels = S(model))
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
    a <- factor(action, levels = seq_along(A(model)), labels = A(model))
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

### check if a field/action has a matrix
.action_is_matrix <- function(model, field, action) {
  f <- model[[field]]
  
  if (is.null(f))
    stop("field ", field, " does not exist in the model!")
  
  if (is.function(f))
    return(FALSE)
  
  m <- model[[field]][[action]]
  if (is.matrix(m) || inherits(m, "Matrix"))
    return(TRUE)
  
  return(FALSE)
}
