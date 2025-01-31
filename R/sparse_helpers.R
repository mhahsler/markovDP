# Internal helpers for sparse representation

# make a matrix sparse if it has low density
.sparsify <- function(x, sparse = TRUE, names = NULL) {
  # we also keep special keywords
  if (is.character(x))
    return(x)
  
  # NULL means as is
  if (is.null(sparse))
    return(x)
  
  if (is.logical(sparse)) {
    if (is.null(dimnames(x)) || any(sapply(dimnames(x), is.null)))
      dimnames(x) <- names
    
    if (!sparse)
      return(as.matrix(x))
    else 
      return(as.csr.matrix(x))
  }
  
  if (is.character(sparse))
    sparse <- match.arg(sparse, c("sparse_no_labels"))
 
  x <- as.csr.matrix(x)
  dimnames(x) <- NULL
  x
}

# accepts a sparse or dense vector or a matrix with one row or column
.sparsify_vector <- function(x, sparse = TRUE, names = NULL) {
  if (is.character(sparse))
    sparse <- match.arg(sparse, c("states", "index", "sparse_no_labels"))
  
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
 
  switch(sparse, 
         sparse_no_labels = as.sparse.vector(x),
         
         states =  {
           if (is(x, "sparseVector")) {
             if (is.null(names))
               stop("state names needed to return states.")
             # faster
             #names[which(x > 0)]
             names[x@i]
           } else
             names(x)[x > 0]
         }, 
         
         index = unname(which(x > 0)), 
        
  )
}


# Translate a probability distribution like the start distribution in MDPs
# Input options:
# * dense/sparse vector
# * dense matrix
# * keyword "uniform"
# * state names
# * state ids
# "-" + state ids to remove
# * neg (removed) ids
#
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
