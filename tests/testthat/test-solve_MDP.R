## context("solve_MDP")

#data(Maze)
#models <- list(Maze)

verbose <- FALSE
#verbose <- TRUE

methods_DP <- c("value_iteration", "policy_iteration", "prioritized_sweeping")
methods_LP <- c("lp")
methods_TD <- c("sarsa", "q_learning", "expected_sarsa")
methods_MC <- c("MC_exploring_starts", "MC_on_policy", "MC_off_policy")
methods_sampling <- c("q_planning")

methods <- c(
  methods_DP,
  # methods_LP, ### to slow
  methods_TD
)

num_states <- length(Maze_orig$states)

times <- data.frame(model = character(0), method = character(0), 
                    user = numeric(0), system = numeric(0),
                    elapsed = numeric(0))


for (model in models_solve) {
  for (m in methods) {
    if (verbose)
      cat("Solving w/", m, ":", model$name, "\n")
    
    t <- system.time(sol <- solve_MDP(model, method = m))
    
    if (verbose)
      cat("time: ", t[3], " sec.\n\n")
    
    times <- rbind(times, data.frame(model = model$name, method = m, 
                                     user = t[1], system = t[2],
                                     elapsed = t[3]))
    
    pol <- policy(sol)
    expect_identical(dim(pol), c(num_states, 3L))
    
    # check_and_fix_MDP(sol)
  }
}

### these methods are slow and need restrictions

for (model in models_solve) {
  for (m in methods_sampling) {
    if (verbose)
      cat("Solving w/", m, ":", model$name,"\n")
    
    t <- system.time(sol <- solve_MDP(model, method = m, n = 10))
    
    if (verbose)
      cat("time: ", t[3], " sec.\n\n")
    times <- rbind(times, data.frame(model = model$name, method = m, 
                                     user = t[1], system = t[2],
                                     elapsed = t[3])) 
    
    pol <- policy(sol)
    expect_identical(dim(pol), c(num_states, 3L))
    
    # check_and_fix_MDP(sol)
  }
}

for (model in models_solve) {
  for (m in methods_MC) {
    if (verbose)
      cat("Solving w/", m, ":", model$name,"\n")
    
    t <- system.time(sol <- solve_MDP(model, method = m, n = 10, horizon = 100))
    
    if (verbose)
      cat("time: ", t[3], " sec.\n\n")
    times <- rbind(times, data.frame(model = model$name, method = m, 
                                     user = t[1], system = t[2],
                                     elapsed = t[3]))
    
    
    pol <- policy(sol)
    expect_identical(dim(pol), c(num_states, 3L))
    
    # check_and_fix_MDP(sol)
  }
}

times[order(times$elapsed), ]
