## simulate_MDP

data(Maze)
verb <- FALSE

# unsolved MDP
simulate_MDP(
  Maze,
  n = 10,
  horizon = 10,
  verbose = verb,
  engine = "r"
)

simulate_MDP(
  Maze,
  n = 10,
  horizon = 10,
  return_states = TRUE,
  verbose = verb,
  engine = "r"
)

M_norm <- normalize_MDP(Maze)
simulate_MDP(
  M_norm,
  n = 10,
  horizon = 10,
  verbose = verb,
  engine = "cpp"
)

simulate_MDP(
  M_norm,
  n = 10,
  horizon = 10,
  return_states = TRUE,
  verbose = verb,
  engine = "cpp"
)

# microbenchmark::microbenchmark(simulate_MDP(Maze, n = 100, horizon = 10, verbose = FALSE, engine = "r"))
# microbenchmark::microbenchmark(simulate_MDP(Maze, n = 100, horizon = 10, verbose = FALSE, engine = "cpp"))

# solved MDP
sol <- solve_MDP(M_norm, discount = 1)

simulate_MDP(sol,
  n = 10,
  horizon = 10,
  verbose = verb
)
