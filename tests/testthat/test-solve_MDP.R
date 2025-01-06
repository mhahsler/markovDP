## context("solve_MDP")

verbose <- interactive()

data(Maze)

#models_solve <- list(Maze)
#models_solve <- list(normalize_MDP(Maze, transition_prob = TRUE, reward = FALSE, sparse = TRUE))
models_solve <- list(normalize_MDP(Maze, transition_prob = TRUE, reward = TRUE, sparse = FALSE))

# precompute dense matrices?
#matrix <- FALSE
matrix <- TRUE

#models_solve <- list(gw_random_maze(20))
#models_solve <- list(normalize_MDP(gw_random_maze(20), transition_prob = TRUE, reward = FALSE, sparse = TRUE))

# gw_plot(models_solve[[1]])


methods_DP <- c("value_iteration", "policy_iteration")
methods_PS <- c("prioritized_sweeping")
methods_LP <- c("lp")
methods_TD <- c("sarsa", "q_learning", "expected_sarsa")
methods_MC <- c("MC_exploring_starts", "MC_on_policy", "MC_off_policy")
methods_sampling <- c("q_planning")

num_states <- length(models_solve[[1]]$states)

# benchmark
#bench <- solve_MDP(models_solve[[1]], method = "lp", discount = 1)
bench <- solve_MDP(models_solve[[1]], method = "value")


timing <- data.frame(model = character(0), 
                     method = character(0),
                     action_discrepancy = numeric(0),
                     time = numeric(0))
solutions <- list()


# need no parameters
for (model in models_solve) {
  for (m in methods_DP) {
    if (verbose)
      cat("Solving w/", m, ":", model$name, "\n")
    
    t <- system.time(sol <- solve_MDP(model, method = m, matrix = matrix))
    
    timing <- rbind(timing, data.frame(model = model$name, 
                                       method = m, 
                                       action_discrepancy = action_discrepancy(sol, bench),
                                       weighted_RMSVE = value_error(sol, bench, weighted = TRUE),
                                       time = t[3]))
    solutions <- append(solutions, setNames(list(sol), m))
    
    if (verbose)
      cat("time: ", t[3], " sec.\n\n")
    
    pol <- policy(sol)
    expect_identical(dim(pol), c(num_states, 3L))
    expect_equal(action_discrepancy(sol, bench), 0)
    
    # check_and_fix_MDP(sol)
  }
}

for (model in models_solve) {
  for (m in methods_PS) {
    for (H_update in c("PS_random", "PS_error", "GenPS")) {
    if (verbose)
      cat("Solving w/", m, ":", model$name, "\n")
    
    t <- system.time(sol <- solve_MDP(model, method = m, H_update = H_update, matrix = matrix))
    
    timing <- rbind(timing, data.frame(model = model$name, 
                                       method = paste0(m, "_", H_update), 
                                       action_discrepancy = action_discrepancy(sol, bench),
                                       weighted_RMSVE = value_error(sol, bench, weighted = TRUE),
                                       time = t[3]))
    solutions <- append(solutions, setNames(list(sol), m))
    
    if (verbose)
      cat("time: ", t[3], " sec.\n\n")
    
    pol <- policy(sol)
    expect_identical(dim(pol), c(num_states, 3L))
    expect_equal(action_discrepancy(sol, bench), 0)
    
    # check_and_fix_MDP(sol)
    }
  }
}


# need no parameters
for (model in models_solve) {
  for (m in methods_LP) {
    if (verbose)
      cat("Solving w/", m, ":", model$name, "\n")
    
    # TD warn about inf horizon
    # has no matrix
    t <- system.time(sol <- solve_MDP(model, method = m, discount = 0.999))
    timing <- rbind(timing, data.frame(model = model$name, 
                                       method = m,
                                       action_discrepancy = action_discrepancy(sol, bench),
                                       weighted_RMSVE = value_error(sol, bench, weighted = TRUE),
                                       time = t[3]))
    solutions <- append(solutions, setNames(list(sol), m))
    
    if (verbose)
      cat("time: ", t[3], " sec.\n\n")
    
    pol <- policy(sol)
    expect_identical(dim(pol), c(num_states, 3L))
    expect_equal(action_discrepancy(sol, bench), 0)
    
    # check_and_fix_MDP(sol)
  }
}



### these methods are slow and need restrictions and may not converge

for (model in models_solve) {
  for (m in methods_sampling) {
    if (verbose)
      cat("Solving w/", m, ":", model$name,"\n")
    
    t <- system.time(sol <- solve_MDP(model, method = m, n = 1000, alpha = 0.1, matrix = matrix))
    
    if (verbose)
      cat("time: ", t[3], " sec.\n\n")
    timing <- rbind(timing, data.frame(model = model$name, 
                                       method = m, 
                                       action_discrepancy = action_discrepancy(sol, bench),
                                       weighted_RMSVE = value_error(sol, bench, weighted = TRUE),
                                       time = t[3]))
    solutions <- append(solutions, setNames(list(sol), m))
     
    pol <- policy(sol)
    expect_identical(dim(pol), c(num_states, 3L))
    
    # check_and_fix_MDP(sol)
  }
}

for (model in models_solve) {
  for (m in c(methods_TD)) {
    if (verbose)
      cat("Solving w/", m, ":", model$name,"\n")
    
    t <- system.time(sol <- solve_MDP(model, method = m, n = 1000, horizon = 100, matrix = matrix))
    
    if (verbose)
      cat("time: ", t[3], " sec.\n\n")
    timing <- rbind(timing, data.frame(model = model$name, 
                                       method = m, 
                                       action_discrepancy = action_discrepancy(sol, bench),
                                       weighted_RMSVE = value_error(sol, bench, weighted = TRUE),
                                       time = t[3]))
    solutions <- append(solutions, setNames(list(sol), m))
    
    pol <- policy(sol)
    expect_identical(dim(pol), c(num_states, 3L))
    
    # check_and_fix_MDP(sol)
  }
}

for (model in models_solve) {
  for (m in c(methods_MC)) {
    if (verbose)
      cat("Solving w/", m, ":", model$name,"\n")
    
    t <- system.time(sol <- solve_MDP(model, method = m, n = 1000, horizon = 100, matrix = matrix))
    
    if (verbose)
      cat("time: ", t[3], " sec.\n\n")
    timing <- rbind(timing, data.frame(model = model$name, 
                                       method = m, 
                                       action_discrepancy = action_discrepancy(sol, bench),
                                       weighted_RMSVE = value_error(sol, bench, weighted = TRUE),
                                       time = t[3]))
    solutions <- append(solutions, setNames(list(sol), m))
    
    pol <- policy(sol)
    expect_identical(dim(pol), c(num_states, 3L))
    
    # check_and_fix_MDP(sol)
  }
}


rownames(timing) <- NULL

#  if (verbose) {
#   library(tidyverse)
#   ggplot(timing, aes(reorder(method, time),
#                      y = time, fill = abbreviate(model))) +
#     geom_bar(stat = "identity")
# }

timing[order(timing$time),]
timing[order(timing$action_discrepancy),]


#cbind(policy(solutions$q_planning), bench = policy(bench)$action)

