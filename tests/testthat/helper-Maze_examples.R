#data("Maze", package = "markovDP")

# defines several versions of the Maze problem
gw <- gw_init(dim = c(3, 4),  
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
T <- function(model, action, start.state, end.state) {
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

R <- rbind(
  R_(                         value = -0.04),
  R_(end.state = "s(2,4)",    value = -1),
  R_(end.state = "s(1,4)",    value = +1),
  R_(start.state = "s(2,4)",  value = 0),
  R_(start.state = "s(1,4)",  value = 0)
)

Maze_orig <- MDP(
  name = "Maze",
  discount = 1,
  horizon = Inf,
  states = gw$states,
  actions = gw$actions,
  transition_prob = T,
  reward = R,
  start = "s(3,1)",
  info = gw$info
)

Maze_orig <- normalize_MDP(
  Maze_orig,
  transition_prob = FALSE,
  reward = FALSE,
  precompute_absorbing = TRUE
)

Maze_dense <- normalize_MDP(Maze_orig, sparse = FALSE)
Maze_sparse <- normalize_MDP(Maze_orig, sparse = TRUE)

Maze_function2 <- Maze_orig
Maze_function2$transition_prob <- function(model, action, start.state) {
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

Maze_function2_sparse <- Maze_function2
Maze_function2_sparse$transition_prob <- function(model, action, start.state) {
  .sparsify_vector(Maze_function2$transition_prob(model, action, start.state))
}

Maze_function2_named <- Maze_function2
Maze_function2_named$transition_prob <- function(model, action, start.state) {
  v <- Maze_function2$transition_prob(model, action, start.state)
  v[v>0]
}

Maze_function3 <- Maze_orig

models <- list(Maze_dense,
#               Maze_sparse,  # sparse has now a different reward where P == 0 is zeroed out
               Maze_function2,
               Maze_function2_sparse,
               Maze_function2_named,
               Maze_function3)


# we test only one since the accessors are tested separately
models_solve <- list(Maze_dense)

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

models <- name_models(models)
#models_solve <- name_models(models_solve)