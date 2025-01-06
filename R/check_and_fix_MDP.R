# Internal function to make sure the definition is complete and
# everything is in the right order and the right factors.

check_and_fix_MDP <- function(x) {
  check_func <- function(x, f, func) {
    # Options are:
    # * 3 formals: "action", "start.state" , "end.state" - returns a single value
    # * 2 formals: "action", "start.state" - returns a vector
    
    req_formals <- c("model", names(formals(func)))
    actual_formals <- names(formals(f))
    
    if (length(actual_formals) == 3L &&
        identical(actual_formals, req_formals[1:3]) ||
        length(actual_formals) == 4L &&
        identical(actual_formals, req_formals[1:4])) {
      return(TRUE)
    }
    
    stop(
      deparse(substitute(f)),
      " has arugments ",
      paste(sQuote(actual_formals), collapse = ", "),
      "\nfunction needs formal arguments: ",
      paste(sQuote(req_formals[1:3]), collapse = ", "),
      " or ",
      paste(sQuote(req_formals[1:4]), collapse = ", ")
    )
  }
  
  check_and_fix_df <- function(x, field, func) {
    # the data.frame needs to match the arguments of the supplied function
    req_columns <- names(formals(func))
    if (is.null(colnames(field))) {
      colnames(field) <- req_columns
    }
    
    if (!identical(colnames(field), req_columns)) {
      stop(
        "The ",
        deparse(substitute(field)),
        " data.frame needs columns named: ",
        paste(sQuote(req_columns), collapse = ", ")
      )
    }
    
    # convert * to NA
    field[field == "*"] <- NA
    field <- type.convert(field, as.is = TRUE)
    
    for (i in grep("action", colnames(field))) {
      if (is.numeric(field[[i]])) {
        field[[i]] <- x$actions[field[[i]]]
      }
      field[[i]] <- as.character(field[[i]])
      
      if (!all(field[[i]] %in% c(x$actions, NA)))
        stop("Unknown action: ", field[[i]][!(field[[i]] %in% c(x$actions, NA))])
      #field[[i]] <- factor(field[[i]], levels = x$actions)
    }
    
    for (i in grep("state", colnames(field))) {
      if (is.numeric(field[[i]])) {
        field[[i]] <- x$states[field[[i]]]
      }
      field[[i]] <- as.character(field[[i]])
      
      if (!all(field[[i]] %in% c(x$states, NA)))
        stop("Unknown action: ", field[[i]][!(field[[i]] %in% c(x$states, NA))])
      #field[[i]] <- factor(field[[i]], levels = x$states)
    }
    # for (i in grep("observation", colnames(field))) {
    #   if (is.numeric(field[[i]])) {
    #     field[[i]] <- x$observations[field[[i]]]
    #   }
    #   field[[i]] <- factor(field[[i]], levels = x$observations)
    # }
    
    field
  }
  
  check_and_fix_list <- function(x, field, func, check_sum = TRUE) {
    if (is.null(names(field))) {
      names(field) <- x$actions
    }
    
    if (all(names(field) != x$actions)) {
      field <- field[x$actions]
    }
    
    for (a in x$actions) {
      if (is.null(field[[a]])) {
        stop(deparse(substitute(field)), ": action ", a, " is missing!")
      }
      
      # we can have a dense or a sparse matrix
      if (is.matrix(field[[a]]) || inherits(x, "sparseMatrix")) {
        if (!identical(dim(field[[a]]), c(length(x$states), length(x$states)))) {
          stop(deparse(substitute(field)), ": matrix for action ",
               a,
               ": has not the right dimensions!")
        }
        if (check_sum && !sum1(field[[a]])) {
          stop(deparse(substitute(field)), ": matrix for action ",
               a,
               ": rows do not add up to 1!")
        }
        if (is.null(dimnames(field[[a]]))) {
          dimnames(field[[a]]) <-
            list(x$states, x$states)
        } else {
          field[[a]][x$states, x$states]
        }
      }
    }
    
    field
  }
  
  ### do the checking
  
  ## expand states, actions and observations if only the number is given
  if (is.numeric(x$states) &&
      length(x$states) == 1L) {
    x$states <- paste0("s", seq_len(x$states))
  }
  
  if (is.numeric(x$actions) &&
      length(x$actions) == 1L) {
    x$actions <- paste0("a", seq_len(x$actions))
  }
  
  ## discount and horizon
  x$discount <- as.numeric(x$discount)
  if (length(x$discount) != 1L ||
      x$discount <= 0 || x$discount > 1) {
    stop("discount has to be a single value in the range (0,1].")
  }
  
  if (is.null(x$horizon)) {
    x$horizon <- Inf
  }
  x$horizon <- as.numeric(x$horizon)
  if (any(x$horizon != floor(x$horizon))) {
    stop("horizon needs to be an integer.")
  }
  
  ## start
  if (is.numeric(x$start) &&
      length(x$start) == length(x$states)) {
    if (!sum1(x$start)) {
      stop("The start probability vector does not add up to 1.")
    }
    if (is.null(names(x$start))) {
      names(x$start) <- x$states
    } else {
      x$start <- x$start[x$states]
    }
  }
  if (any(is.na(x$start))) {
    stop("start containes undefined start states.")
  }
  if (is.character(x$start)) {
    if (!(identical(x$start, "uniform") ||
          all(x$start %in% x$states))) {
      stop(
        "when using characters for start, then it needs to be the keyword 'uniform' or a set of start states."
      )
    }
  }
  
  ## transitions
  if (is.null(x$transition_prob)) {
    stop("transition_prob cannot be missing!")
  }
  
  if (is.function(x$transition_prob)) {
    check_func(x, x$transition_prob, P_)
  }
  
  else if (is.data.frame(x$transition_prob)) {
    x$transition_prob <- check_and_fix_df(x, x$transition_prob, P_)
  }
  
  # list of matrices or keywords
  else {
    x$transition_prob <- check_and_fix_list(x, x$transition_prob, P_, 
                                            check_sum = TRUE)
  }
  
  ## reward
  if (is.null(x$reward)) {
    stop("reward cannot be missing!")
  }
  
  if (is.function(x$reward)) {
    check_func(x, x$reward, R_)
  }
  
  else if (is.data.frame(x$reward)) {
    x$reward <- check_and_fix_df(x, x$reward, R_)
  }
  
  # list of state x state matrices
  else {
    x$reward <- check_and_fix_list(x, x$reward, R_, check_sum = FALSE)
  }
  
  
  ## MDP has no terminal values
  if (!is.null(x$terminal_values)) {
    stop("MDPs do not have terminal_values!")
  }
  
  ### check solution
  if (!is.null(x$solution)) {
    if (is.null(policy)) {
      stop("Policy in solved MDP missing!")
    }
    x$solution$policy <- lapply(
      x$solution$policy,
      FUN = function(p) {
        if (is.null(p$action) ||
            length(p$action) != length(x$states) ||
            !all(p$action %in% x$actions)) {
          stop("Malformed action in MDP policy solution.")
        }
        if (is.null(p$state)) {
          p <- cbind(state = x$states, p)
        }
        p
      }
    )
  }
  
  x
}
