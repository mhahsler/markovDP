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
    value_matrix(
      x,
      "reward",
      action,
      start.state,
      end.state,
      sparse,
      FALSE
    )
  }


# TODO: This is currently unused!
# try to convert a reward matrix into a sparse data.frame
reward_list2df <- function(x) {
  message("This is experimental!")
  actions <- names(x)
  states <- names(x[[1L]])
  
  reward_df <- data.frame()
  for (action in seq_along(actions)) {
    for (start.state in seq_along(states)) {
      for (end.state in seq_along(states)) {
        reward_df <- rbind(reward_df,
                           R_(action, start.state, end.state, value = x[[action]][start.state, end.state]))
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
  
  rownames(reward_df) <- NULL
  
  reward_df
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
