data("Maze", package = "markovDP")

# defines several versions of the Maze problem
gw <- gw_init(
  dim = c(3, 4),
  blocked_states = "s(2,2)",
  absorbing_states = c("s(1,4)", "s(2,4)"),
  state_labels = list(
    "s(3,1)" = "Start",
    "s(2,4)" = "-1",
    "s(1,4)" = "Goal: +1"
  )
)

# the transition function is stochastic so we cannot use the standard
# gridworld gw$transition_prob() function and have to replace it
P <- function(model, action, start.state, end.state) {
  action <- match.arg(action, choices = model$actions)
  
  # absorbing states
  if (start.state %in% model$info$absorbing_states) {
    if (start.state == end.state) {
      return(1)
    } else {
      return(0)
    }
  }
  
  if (action %in% c("up", "down")) {
    error_direction <- c("right", "left")
  } else {
    error_direction <- c("up", "down")
  }
  
  rc <- gw_s2rc(start.state)
  delta <- list(
    up = c(-1, 0),
    down = c(+1, 0),
    right = c(0, +1),
    left = c(0, -1)
  )
  P <- matrix(0, nrow = 3, ncol = 4)
  
  add_prob <- function(P, rc, a, value) {
    new_rc <- rc + delta[[a]]
    # stay in place for moved to a non-existing state
    if (!(gw_rc2s(new_rc) %in% model$states)) {
      new_rc <- rc
    }
    P[new_rc[1], new_rc[2]] <- P[new_rc[1], new_rc[2]] + value
    P
  }
  
  P <- add_prob(P, rc, action, .8)
  P <- add_prob(P, rc, error_direction[1], .1)
  P <- add_prob(P, rc, error_direction[2], .1)
  
  P[rbind(gw_s2rc(end.state))]
}

# Transitions with 2 parameters (is now standard)
P2 <- function(model, action, start.state) {
  action <- match.arg(action, choices = model$actions)
  
  P <- structure(numeric(length(model$states)), names = model$states)
  
  # absorbing states
  if (start.state %in% model$info$absorbing_states) {
    P[start.state] <- 1
    return(P)
  }
  
  if (action %in% c("up", "down")) {
    error_direction <- c("right", "left")
  } else {
    error_direction <- c("up", "down")
  }
  
  rc <- gw_s2rc(start.state)
  delta <- list(
    up = c(-1, 0),
    down = c(+1, 0),
    right = c(0, +1),
    left = c(0, -1)
  )
  
  # there are 3 directions. For blocked directions, stay in place
  # 1) action works .8
  rc_new <- gw_rc2s(rc + delta[[action]])
  if (rc_new %in% model$states)
    P[rc_new] <- .8
  else
    P[start.state] <- .8
  
  # 2) off to the right .1
  rc_new <- gw_rc2s(rc + delta[[error_direction[1]]])
  if (rc_new %in% model$states)
    P[rc_new] <- .1
  else
    P[start.state] <-  P[start.state] + .1
  
  # 3) off to the left .1
  rc_new <- gw_rc2s(rc + delta[[error_direction[2]]])
  if (rc_new %in% model$states)
    P[rc_new] <- .1
  else
    P[start.state] <-  P[start.state] + .1
  
  P
}

R <- rbind(
  R_(value = -0.04),
  R_(end.state = "s(2,4)", value = -1),
  R_(end.state = "s(1,4)", value = +1),
  R_(start.state = "s(2,4)", value = 0),
  R_(start.state = "s(1,4)", value = 0)
)

R_func <- function(model, action, start.state, end.state) {
  if (start.state %in% c("s(2,4)", "s(1,4)"))
    return(0)
  
  if (end.state == "s(2,4)")
    return(-1)
  if (end.state == "s(1,4)")
    return(+1)
  
  return(-0.04)
}

Maze_function3 <- MDP(
  name = "Maze",
  discount = 1,
  horizon = Inf,
  states = gw$states,
  actions = gw$actions,
  transition_prob = P,
  reward = R,
  start = "s(3,1)",
  info = gw$info
)

Maze_function3 <- normalize_MDP(
  Maze_function3,
  transition_prob = FALSE,
  reward = FALSE,
  precompute_absorbing = TRUE
)

Maze_dense <- normalize_MDP(Maze_function3, sparse = FALSE)
Maze_sparse <- normalize_MDP(Maze_function3, sparse = TRUE)

Maze_function2 <- Maze_function3
Maze_function2$transition_prob <- P2

# function returns a sparse vector
Maze_function2_sparse <- Maze_function2
Maze_function2_sparse$transition_prob <- function(model, action, start.state) {
  .sparsify_vector(Maze_function2$transition_prob(model, action, start.state))
}

# function returns a named vector
Maze_function2_named <- Maze_function2
Maze_function2_named$transition_prob <- function(model, action, start.state) {
  v <- Maze_function2$transition_prob(model, action, start.state)
  v[v > 0]
}

# original has a dens transition matrix and a data frame for rewards
Maze_orig <- Maze_function2
Maze_orig <- normalize_MDP(Maze_function2, transition_prob = TRUE, reward = FALSE)

Maze_reward_function <- Maze_orig
Maze_reward_function$reward <- R_func
Maze_reward_trans_function <- Maze_reward_function
Maze_reward_trans_function$transition_prob <- P2

# test lists
models_matrix <- list(
  dense = Maze_dense, 
  sparse = Maze_sparse
  )

models_trans_function <- list(
  f2 = Maze_function2,
  f2_sparse = Maze_function2_sparse,
  f2_named = Maze_function2_named,
  f3 = Maze_function3
  )

models_reward_function <- list(
  f = Maze_reward_function,
  ff = Maze_reward_trans_function
  )

name_models <- function(models)
  lapply(
    models,
    FUN = function(m) {
      m$name <- paste(
        m$name,
        "- T:",
        class(m$transition_prob),
        if (is.list(m$transition_prob))
          paste0("(", paste(
            sapply(
              m$transition_prob,
              FUN = function(x)
                class(x)[1]
            ), collapse = ", "
          ), ")"),
        if (is.function(m$transition_prob))
          paste0("(", paste(names(
            formals(m$transition_prob)
          ), collapse = ", "), ")"),
        "- R:",
        class(m$reward),
        if (is.list(m$reward))
          paste0("(", paste(
            sapply(
              m$reward,
              FUN = function(x)
                class(x)[1]
            ), collapse = ", "
          ), ")"),
        "- SOLVED:",
        is_solved_MDP(m)
      )
      m
    }
  )

models_matrix <- name_models(models_matrix)
models_trans_function <- name_models(models_trans_function)