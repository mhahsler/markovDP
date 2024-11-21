## context("solve_MDP")

# models_solve <- list(gw_random_maze(20))
# gw_plot(models_solve[[1]])

verbose <- interactive()

methods_DP <- c("value_iteration", "policy_iteration")
methods_PS <- c("prioritized_sweeping")
methods_LP <- c("lp")
methods_TD <- c("sarsa", "q_learning", "expected_sarsa")
methods_MC <- c("MC_exploring_starts", "MC_on_policy", "MC_off_policy")
methods_sampling <- c("q_planning")

num_states <- length(models_solve[[1]]$states)
bench <- solve_MDP(models_solve[[1]], method = "lp", discount = 0.999)


timing <- data.frame(model = character(0), 
                     method = character(0),
                     action_discrepancy = numeric(0),
                     time = numeric(0))

solutions <- list()

# benchmark
bench <- solve_MDP(models_solve[[1]], method = "value")

# need no parameters
for (model in models_solve) {
  for (m in methods_DP) {
    if (verbose)
      cat("Solving w/", m, ":", model$name, "\n")
    
    t <- system.time(sol <- solve_MDP(model, method = m))
    
    timing <- rbind(timing, data.frame(model = model$name, 
                                       method = m, 
                                       action_discrepancy = action_discrepancy(sol, bench),
                                       time = t[3]))
    solutions <- append(solutions, list(sol))
    
    if (verbose)
      cat("time: ", t[3], " sec.\n\n")
    
    pol <- policy(sol)
    expect_identical(dim(pol), c(num_states, 3L))
    
    # check_and_fix_MDP(sol)
  }
}

for (model in models_solve) {
  for (m in methods_PS) {
    for (H_update in c("PS_random", "PS_error", "GenPS")) {
    if (verbose)
      cat("Solving w/", m, ":", model$name, "\n")
    
    t <- system.time(sol <- solve_MDP(model, method = m, H_update = H_update))
    
    timing <- rbind(timing, data.frame(model = model$name, 
                                       method = paste0(m, "_", H_update), 
                                       action_discrepancy = action_discrepancy(sol, bench),
                                       time = t[3]))
    solutions <- append(solutions, list(sol))
    
    if (verbose)
      cat("time: ", t[3], " sec.\n\n")
    
    pol <- policy(sol)
    expect_identical(dim(pol), c(num_states, 3L))
    
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
    t <- system.time(suppressWarnings((sol <- solve_MDP(model, method = m))))
    timing <- rbind(timing, data.frame(model = model$name, 
                                       method = m,
                                       action_discrepancy = action_discrepancy(sol, bench),
                                       time = t[3]))
    solutions <- append(solutions, list(sol))
    
    if (verbose)
      cat("time: ", t[3], " sec.\n\n")
    
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
    
    t <- system.time(sol <- solve_MDP(model, method = m, n = 10000))
    
    if (verbose)
      cat("time: ", t[3], " sec.\n\n")
    timing <- rbind(timing, data.frame(model = model$name, 
                                       method = m, 
                                       action_discrepancy = action_discrepancy(sol, bench),
                                       time = t[3]))
    solutions <- append(solutions, list(sol))
   
   
     
    pol <- policy(sol)
    expect_identical(dim(pol), c(num_states, 3L))
    
    # check_and_fix_MDP(sol)
  }
}

for (model in models_solve) {
  for (m in c(methods_TD)) {
    if (verbose)
      cat("Solving w/", m, ":", model$name,"\n")
    
    t <- system.time(sol <- solve_MDP(model, method = m, n = 1000, horizon = 100))
    
    if (verbose)
      cat("time: ", t[3], " sec.\n\n")
    timing <- rbind(timing, data.frame(model = model$name, 
                                       method = m, 
                                       action_discrepancy = action_discrepancy(sol, bench),
                                       time = t[3]))
    
    
    solutions <- append(solutions, list(sol))
    
    pol <- policy(sol)
    expect_identical(dim(pol), c(num_states, 3L))
    
    # check_and_fix_MDP(sol)
  }
}

for (model in models_solve) {
  for (m in c(methods_MC)) {
    if (verbose)
      cat("Solving w/", m, ":", model$name,"\n")
    
    t <- system.time(sol <- solve_MDP(model, method = m, n = 1000, horizon = 100))
    
    if (verbose)
      cat("time: ", t[3], " sec.\n\n")
    timing <- rbind(timing, data.frame(model = model$name, 
                                       method = m, 
                                       action_discrepancy = action_discrepancy(sol, bench),
                                       time = t[3]))
    
    
    solutions <- append(solutions, list(sol))
    
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
