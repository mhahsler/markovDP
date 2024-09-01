## context("MDP_functions")

data("Maze")
Maze_dense <- normalize_MDP(Maze, sparse = FALSE)
Maze_sparse <- normalize_MDP(Maze, sparse = TRUE)
Maze_function <- Maze

Maze_function$transition_prob <- function(action, start.state, end.state) {
  actions <- c("up", "right", "down", "left")
  states <- c(
    "s(1,1)", "s(2,1)", "s(3,1)", "s(1,2)", "s(3,2)", "s(1,3)",
    "s(2,3)", "s(3,3)", "s(1,4)", "s(2,4)", "s(3,4)"
  )
  
  action <- match.arg(action, choices = actions)
  
  # absorbing states
  if (start.state %in% c("s(1,4)", "s(2,4)")) {
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
  
  rc <- gridworld_s2rc(start.state)
  delta <- list(
    up = c(-1, 0),
    down = c(+1, 0),
    right = c(0, +1),
    left = c(0, -1)
  )
  P <- matrix(0, nrow = 3, ncol = 4)
  
  add_prob <- function(P, rc, a, value) {
    new_rc <- rc + delta[[a]]
    if (!(gridworld_rc2s(new_rc) %in% states)) {
      new_rc <- rc
    }
    P[new_rc[1], new_rc[2]] <- P[new_rc[1], new_rc[2]] + value
    P
  }
  
  P <- add_prob(P, rc, action, .8)
  P <- add_prob(P, rc, error_direction[1], .1)
  P <- add_prob(P, rc, error_direction[2], .1)
  P[rbind(gridworld_s2rc(end.state))]
}

# transitions
tr <- transition_matrix(Maze, "up", 1)
tr_dense <- transition_matrix(Maze_dense, "up", 1)
tr_sparse <- transition_matrix(Maze_sparse, "up", 1)
tr_function <- transition_matrix(Maze_function, "up", 1)

expect_equal(as.vector(tr), as.vector(tr_dense))
expect_equal(tr, tr_sparse)
expect_equal(tr_function, tr_dense)

tr <- transition_matrix(Maze, sparse = FALSE)
tr_dense <- transition_matrix(Maze_dense, sparse = FALSE)
tr_sparse <- transition_matrix(Maze_sparse, sparse = FALSE)
tr_function <- transition_matrix(Maze_function, sparse = FALSE)

expect_equal(tr, tr_dense)
expect_equal(tr, tr_sparse)
expect_equal(tr, tr_function)

tr <- transition_matrix(Maze, sparse = TRUE)
tr_dense <- transition_matrix(Maze_dense, sparse = TRUE)
tr_sparse <- transition_matrix(Maze_sparse, sparse = TRUE)
tr_function <- transition_matrix(Maze_function, sparse = TRUE)

expect_equal(tr, tr_dense)
expect_equal(tr, tr_sparse)
expect_equal(tr, tr_function)

# reward
v <- reward_matrix(Maze, "up", 1, 1)
v_dense <- reward_matrix(Maze_dense, "up", 1, 1)
v_sparse <- reward_matrix(Maze_sparse, "up", 1, 1)

expect_equal(v, v_dense)
expect_equal(v, v_sparse)

# absorbing states
s_abs <- c("s(1,4)", "s(2,4)")
expect_equal(names(which(absorbing_states(Maze))), s_abs)
expect_equal(names(which(absorbing_states(Maze_dense))), s_abs)
expect_equal(names(which(absorbing_states(Maze_sparse))), s_abs)
expect_equal(names(which(absorbing_states(Maze_function))), s_abs)

# reachable states
s_reach <- Maze$states
expect_equal(names(which(reachable_states(Maze))), s_reach)
expect_equal(names(which(reachable_states(Maze_dense))), s_reach)
expect_equal(names(which(reachable_states(Maze_sparse))), s_reach)
expect_equal(names(which(reachable_states(Maze_function))), s_reach)
