## simulate_MDP

verbose <- FALSE
verb <- FALSE

#verbose <- TRUE
#verb <- TRUE

# for unsolved model, sample uses a randomized policy 
#    (i.e., an epsilon-soft policy with epsilon = 1).
# add some solved models where sample used the optimal policy

Maze_sparse_sol <- solve_MDP(Maze_sparse)
Maze_function2_sol <- solve_MDP(Maze_function2)  ### will fallback to R

models_solved <- list(Maze_sparse_sol,
                      Maze_function2_sol)


models_solved <- 
  name_models(models_solved)

for (m in c(models_solve, models_solved)) {
  if (verbose)
    cat("\n\nSampling w/R from: ", m$name, " \n")
  
  n <- 10
  
  ret <- sample_MDP(
    m,
    n = n,
    horizon = 10,
    trajectories = TRUE,
    verbose = verb,
    engine = "r"
  )
  
  expect_length(ret$reward, n)
  expect_length(ret$action_cnt, length(m$actions))
  expect_length(ret$state_cnt, length(m$states))
  expect_equal(ncol(ret$trajectories), 6L)
  
  
  if (verbose)
    cat("Sampling w/cpp from: ", m$name, "\n")
  
  ret <- sample_MDP(
    m,
    n = n,
    horizon = 10,
    verbose = verb,
    trajectories = TRUE,
    engine = "cpp"
  )
  
  expect_length(ret$reward, n)
  expect_length(ret$action_cnt, length(m$actions))
  expect_length(ret$state_cnt, length(m$states))
  expect_equal(ncol(ret$trajectories), 6L)
}


# microbenchmark::microbenchmark(simulate_MDP(Maze, n = 100, horizon = 10, verbose = FALSE, engine = "r"))
# microbenchmark::microbenchmark(simulate_MDP(Maze, n = 100, horizon = 10, verbose = FALSE, engine = "cpp"))

