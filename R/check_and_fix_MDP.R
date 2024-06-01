# Internal function to make sure the definition is complete and
# everything is in the right order and the right factors
check_and_fix_MDP <- function(x) {
  check_func <- function(x, func, name) {
    req_formals <- head(names(formals(func)), -1)
    if (!identical(names(formals(x)), req_formals)) {
      stop(
        name,
        " function needs formal arguments: ",
        paste(sQuote(req_formals), collapse = ", ")
      )
    }
  }

  check_df <- function(x, field, func) {
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
      field[[i]] <- factor(field[[i]], levels = x$actions)
    }

    for (i in grep("state", colnames(field))) {
      if (is.numeric(field[[i]])) {
        field[[i]] <- x$states[field[[i]]]
      }
      field[[i]] <- factor(field[[i]], levels = x$states)
    }
    # for (i in grep("observation", colnames(field))) {
    #   if (is.numeric(field[[i]])) {
    #     field[[i]] <- x$observations[field[[i]]]
    #   }
    #   field[[i]] <- factor(field[[i]], levels = x$observations)
    # }

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

  # if (inherits(x, "POMDP")) {
  #   if (is.numeric(x$observations) &&
  #     length(x$observations) == 1L) {
  #     x$observations <- paste0("o", seq_len(x$observations))
  #   }
  # }

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
    if (!(identical(x$start, "uniform") || all(x$start %in% x$states))) {
      stop(
        "when using characters for start, then it needs to be the keyword 'uniform' or a set of start states."
      )
    }
  }

  ## transitions
  if (is.null(x$transition_prob)) {
    stop("transition_prob cannot be missing!")
  }

  # if we have matrices then check and add names
  if (is.function(x$transition_prob)) {
    check_func(x$transition_prob, T_, "transition_prob")
  } # x$transition_prob <- transition_matrix(x)
  else if (is.data.frame(x$transition_prob)) {
    x$transition_prob <- check_df(x, x$transition_prob, T_)
  } else {
    # action names and order
    if (is.null(names(x$transition_prob))) {
      names(x$transition_prob) <- x$actions
    }
    if (all(names(x$transition_prob) != x$actions)) {
      x$transition_prob <- x$transition_prob[x$actions]
    }

    for (a in x$actions) {
      if (is.null(x$transition_prob[[a]])) {
        stop("transition_prob for action ", a, " is missing!")
      }
      if (is.matrix(x$transition_prob[[a]])) {
        if (!identical(dim(x$transition_prob[[a]]), c(length(x$states), length(x$states)))) {
          stop(
            "transition_prob matrix for action ",
            a,
            ": has not the right dimensions!"
          )
        }
        if (!sum1(x$transition_prob[[a]])) {
          stop(
            "transition_prob matrix for action ",
            a,
            ": rows do not add up to 1!"
          )
        }
        if (is.null(dimnames(x$transition_prob[[a]]))) {
          dimnames(x$transition_prob[[a]]) <-
            list(x$states, x$states)
        } else {
          x$transition_prob[[a]][x$states, x$states]
        }
      }
    }
  }

  ## reward
  if (is.null(x$reward)) {
    stop("reward cannot be missing!")
  }

  # MDP has no observations
  R_ <- function(action = NA,
                 start.state = NA,
                 end.state = NA,
                 value) {
    data.frame(
      action = action,
      start.state = start.state,
      end.state = end.state,
      value = as.numeric(value),
      stringsAsFactors = FALSE
    )
  }

  if (is.function(x$reward)) {
    check_func(x$reward, R_, "reward")
  } else if (is.data.frame(x$reward)) {
    x$reward <- check_df(x, x$reward, R_)
  } else {
    if (is.null(names(x$reward))) {
      names(x$reward) <- x$actions
    }
    if (all(names(x$reward) != x$actions)) {
      x$reward <- x$reward[x$actions]
    }

    for (a in x$actions) {
      if (is.null(x$reward[[a]])) {
        stop("reward for action ", a, " is missing!")
      }
      for (s in x$states) {
        if (is.null(x$reward[[a]][[s]])) {
          stop(
            "reward for action ",
            a,
            " and state ",
            s,
            " is missing!"
          )
        }
        if (is.matrix(x$reward[[a]][[s]])) {
          if (!identical(dim(x$reward[[a]][[s]]), c(length(x$states), length(x$observations)))) {
            stop(
              "reward matrix for action ",
              a,
              " and start.state ",
              s,
              ": does not have the right dimensions!"
            )
          }
          if (is.null(dimnames(x$reward[[a]][[s]]))) {
            dimnames(x$reward[[a]][[s]]) <-
              list(x$states, NULL)
          } else {
            x$reward[[a]][[s]][x$states, NULL]
          }
        }
      }
    }
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
