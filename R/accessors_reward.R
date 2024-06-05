# Accessor Functions for Rewards
#
# Representations:
# Default:
# * Sparse (df): A data.frame with: action start.state end.state observation value
#
# Others:
# * Dense (list): A action list -> start.state list -> end.state x observation matrix
# * A function can be converted to a list


### Note: Sparse is a data.frame for rewards
#' @include accessors.R
#' @rdname accessors
#' @export
reward_matrix <- function(x,
                          action = NULL,
                          start.state = NULL,
                          end.state = NULL,
                          ...,
                          sparse = FALSE) {
  UseMethod("reward_matrix")
}

#' @export
reward_matrix.MDP <-
  function(x,
           action = NULL,
           start.state = NULL,
           end.state = NULL,
           ...,
           state_matrix = TRUE,
           sparse = FALSE) {
    reward <- x[["reward"]]

    if (is.null(action) && (!is.null(start.state) || !is.null(end.state))) {
      stop("action needs to be specified if states are specified!")
    }
    if (!is.null(action) && is.null(start.state) && (!is.null(end.state))) {
      stop("start.state needs to be specified if end.state is specified!")
    }

    # convert functions first
    if (is.function(reward)) {
      # shortcut for a single value
      if (!is.null(action) && !is.null(start.state) && !is.null(end.state)) {
        if (is.numeric(action)) action <- x$actions[action]
        if (is.numeric(start.state)) start.state <- x$states[start.state]
        if (is.numeric(end.state)) end.state <- x$states[end.state]

        return(reward(action, start.state, end.state))
      }

      reward <- reward_function2list(reward, x)
    }

    # return as is
    if (is.null(sparse)) {
      return(reward)
    }

    # data.frame
    if (is.data.frame(reward)) {
      if (sparse) {
        return(reward)
      } else {
        return(reward_df2value(reward, action, start.state, end.state, 
                               state_matrix = state_matrix))
      }
    }

    # we have a list of list of  matrices
    if (sparse) {
      return(reward_list2df(reward))
    }

    # subset
    reward_list2value(reward, action, start.state, end.state)
  }



### NOTE: this does not check if the actions/states/observations are valid names!
reward_df2value <-
  function(df,
           action = NULL,
           start.state = NULL,
           end.state = NULL,
           observation = NULL,
           state_matrix = FALSE, # MDP has no observations
           sparse = FALSE) {
    actions <- levels(df$action)
    states <- levels(df$end.state)

    # MDP has no observations

    observations <- NA
    observation <- 1L

    if (is.null(action)) {
      l <- sapply(
        actions,
        FUN = function(a) {
          sapply(
            states,
            FUN = function(s) {
              .sparsify(reward_df2value(df, a, s), sparse = sparse)
            },
            simplify = FALSE
          )
        },
        simplify = FALSE
      )
      
      if (state_matrix) {
        l <- sapply(actions, FUN = function(a) {
          r <- t(do.call(cbind, l[[a]]))
          rownames(r) <- states
          r
        }, simplify = FALSE)
        
      }
      
      return(l)
    }

    if (is.null(start.state)) {
      l <-
        sapply(
          states,
          FUN = function(s) {
            .sparsify(reward_df2value(df, action, s), sparse = sparse)
          },
          simplify = FALSE
        )

      if (state_matrix) {
        l <- t(do.call(cbind, l))
        rownames(l) <- states
      }
      
      return(l)
    }

    if (is.numeric(action)) {
      action <- actions[action]
    }
    if (is.numeric(start.state)) {
      start.state <- states[start.state]
    }


    if (is.null(end.state)) {
      # end.state vector
      df <- df[
        (is.na(df$action) | df$action == action) &
          (is.na(df$start.state) |
            df$start.state == start.state), ,
        drop = FALSE
      ]

      v <- structure(numeric(length(states)), names = states)

      for (i in seq_len(nrow(df))) {
        e.s <- df$end.state[i]

        if (is.na(e.s)) {
          e.s <- states
        }

        v[e.s] <- df$value[i]
      }

      return(v)
    }

    # value
    if (is.numeric(end.state)) {
      end.state <- states[end.state]
    }

    val <- df$value[(is.na(df$action) | df$action == action) &
      (is.na(df$start.state) |
        df$start.state == start.state) &
      (is.na(df$end.state) |
        df$end.state == end.state)]

    if (length(val) == 0L) {
      return(0)
    }

    return(tail(val, 1L))
  }

reward_df2list <- function(df, sparse = FALSE) {
  reward_df2value(df, sparse = sparse)
}

reward_list2df <- function(x) {
  actions <- names(x)
  states <- names(x[[1L]])
  observations <- colnames(x[[1L]][[1L]])

  reward_df <- data.frame()
  for (action in seq_along(actions)) {
    for (start.state in seq_along(states)) {
      for (observation in seq_along(observations)) {
        for (end.state in seq_along(states)) {
          reward_df <- rbind(
            reward_df,
            R_(
              action,
              start.state,
              end.state,
              observation,
              value = x[[action]][[start.state]][end.state, observation]
            )
          )
        }
      }
    }
  }

  # use most common as the default
  # TODO: more simplifications
  tbl <- table(reward_df$value)
  most_common <- as.numeric(names(tbl[which.max(tbl)]))
  default_lines <- R_(value = most_common)
  reward_df <-
    reward_df[reward_df$value != most_common, , drop = FALSE]

  reward_df <- rbind(default_lines, reward_df)
  ## TODO: MDP

  reward_df[["action"]] <-
    factor(reward_df[["action"]], labels = actions, levels = seq_along(actions))
  reward_df[["start.state"]] <-
    factor(reward_df[["start.state"]], labels = states, levels = seq_along(states))
  reward_df[["end.state"]] <-
    factor(reward_df[["end.state"]], labels = states, levels = seq_along(states))
  reward_df[["observation"]] <-
    factor(reward_df[["observation"]], labels = observations, levels = seq_along(observations))


  rownames(reward_df) <- NULL

  reward_df
}


reward_function2list <- function(f, model, sparse = FALSE) {
  actions <- model$actions
  states <- model$states
  observations <- model$observations

  f <- Vectorize(f)

  sapply(
    actions,
    FUN = function(a) {
      sapply(
        states,
        FUN = function(s) {
          p <- outer(
            states,
            observations,
            FUN = function(end, o) {
              f(
                action = a,
                start.state = s,
                end.state = end,
                observation = o
              )
            }
          )
          dimnames(p) <- list(states, observations)
          p <- .sparsify(p, sparse)
        },
        simplify = FALSE
      )
    },
    simplify = FALSE
  )
}


### this just subsets the matrix list
reward_list2value <-
  function(m,
           action = NULL,
           start.state = NULL,
           end.state = NULL) {

    if (is.null(action)) {
      return(m)
    }
    
    if (is.null(start.state)) {
      return(m[[action]])
    }

    if (is.null(end.state)) {
      return(m[[action]][start.state, ])
    }

    return(m[[action]][start.state, end.state])
  }

.reward_range <- function(model) {
  reward <- model[["reward"]]
  
  if (is.data.frame(reward)) 
    return (range(reward$value))
  
  # function is slow!
  if (is.function(reward)) 
    reward <- reward_function2list(reward, model)
 
  # it is a list of matrices 
  range(unlist(reward, recursive = TRUE)) 
}
