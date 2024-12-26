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
#' @rdname accessors
#' @export
reward_matrix <- function(model,
                          action = NULL,
                          start.state = NULL,
                          end.state = NULL,
                          ...,
                          sparse = NULL,
                          simplify = FALSE) {
  UseMethod("reward_matrix")
}

#' @export
reward_matrix.MDP <-
  function(model,
           action = NULL,
           start.state = NULL,
           end.state = NULL,
           ...,
           state_matrix = TRUE,
           sparse = NULL,
           simplify = FALSE) {
    value_matrix(
      model,
      "reward",
      action,
      start.state,
      end.state,
      sparse,
      simplify,
      trans_keyword = FALSE
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

## fast way to get the range from a reward data.frame
.reward_range <- function(model) {
  
  # this is fast
  if (is.data.frame(model$reward))
    return (range(model$reward$value))
  
  # converting a function is slow!
  reward <- reward_matrix(model, sparse = NULL)
  
  c(min(sapply(reward, min)), max(sapply(reward, max)))
}
