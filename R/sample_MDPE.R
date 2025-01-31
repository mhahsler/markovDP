#' Sample Trajectories from an MDPE
#'
#' Sample trajectories through using a MDPE.
#'
#' @family MDPE

#' @importFrom stats runif
#'
#' @param model an MDPE model.
#' @param n number of trajectories.
#' @param start start state.
#' @param horizon epochs end once an absorbing state is reached or after
#'  the maximal number of epochs specified via `horizon`. If `NULL` then the
#'  horizon for the model is used.
#' @param epsilon the probability of random actions for using an epsilon-greedy policy.
#'  Default for solved models is 0 and for unsolved model 1.
#' @param trajectories logical; return the complete trajectories.
#' @param progress show a progress bar?
#' @param verbose report used parameters
#' @param ... further arguments are ignored.
#' @return A list with elements:
#'  * `avg_reward`: The average discounted reward.
#'  * `reward`: Reward for each trajectory.
#'  * `trajectories`: A data.frame with the trajectories. Each row
#'    contains the `episode` id, the `time` step, the state `s`,
#'    the chosen action `a`,
#'    the reward `r`, and the next state `s_prime`. Trajectories are
#'    only returned for `trajectories = TRUE`.
#' @author Michael Hahsler
#' @examples
#' # enable parallel simulation
#' # doParallel::registerDoParallel()
#'
#' # Create a simple maze with the layout:
#' # XXXXXX
#' # XSX  X  
#' # X    X
#' # X  X X
#' # X  XGX
#' # XXXXXX
#' 
#' model <- gw_maze_MDPE(
#'            dim = c(4, 4),
#'            start = c(1, 1),
#'            goal = c(4, 4),
#'            walls = list(c(1, 2), c(3, 3), c(3, 4)),
#'            discount = 0.95,
#'            name = "Simple Maze"
#'        )
#' model
#' 
#' sim <- sample_MDP(model, horizon = 500, n = 1, 
#'                    verbose = TRUE, trajectories = TRUE)
#' sim
#' @export
sample_MDP.MDPE <-
  function(model,
           n,
           start = NULL,
           horizon = NULL,
           epsilon = NULL,
           trajectories = FALSE,
           progress = TRUE,
           verbose = FALSE,
           ...) {
    .nodots(...)
    
    #solved <- is_solved_MDP(model)
    warning("Fix solved with policy!")
    solved <- FALSE
    
    n <- as.integer(n)
    
    start <- start %||% model$start
    
    horizon <- horizon %||% model$horizon %||% Inf
    if (is.infinite(horizon))
      stop("Finite sampling horizon is needed!")
    horizon <- as.integer(horizon)
    
    epsilon <- epsilon %||% ifelse(solved, 0, 1)
    
    if (!solved && epsilon != 1) {
      stop("epsilon has to be 1 for unsolved models.")
    }
    
    disc <- model$discount %||% 1
    
    actions <- A(model)
    
    if (verbose) {
      cat("Sampling MDPE trajectories.\n")
      cat("- horizon:", horizon, "\n")
      cat("- n:",
          n,
          "- parallel workers:",
          foreach::getDoParWorkers(),
          "\n")
      cat("- epsilon:", epsilon, "\n")
      cat("- discount factor:", disc, "\n")
      cat("\n")
    }
    
    # Progressbar does not work with foreach
    if (foreach::getDoParWorkers() != 1L)
      progress <- FALSE
    
    if (progress)
      pb <- my_progress_bar(n, name = "sample_MDPE")
   
    #warning("Debug mode on!!!")
    #sim <- for(i in 1:n){
    sim <- foreach(i = 1:n) %dopar% {
      if (progress)
        pb$tick()
      
      action_cnt <- rep(0L, length(actions))
      names(action_cnt) <- actions
      
      state_cnt <- fastmap()
      rew <- 0
      
      # find a initial state
      s <- start
      
      if (trajectories) {
        trajectory <- data.frame(
          episode = rep(NA_integer_, horizon),
          time = rep(NA_integer_, horizon),
          s = NA_character_,
          a = NA_integer_,
          r = NA_real_,
          s_prime = NA_character_
        )
      } else {
        trajectory <- NULL
      }
      
      for (j in seq_len(horizon)) {
        # epsilon soft policy
        #if (runif(1) < epsilon) {
          a <- sample.int(length(actions), 1L, replace = TRUE)
        #} else {
        #  a <- pol[[.get_pol_index(model, j)]][s]
        #}
        
        action_cnt[a] <- action_cnt[a] + 1L
        s_label <- features2state(s)
        state_cnt$set(s_label, state_cnt$get(s_label, missing = 0L) + 1L)
        
        s_prev <- s
        result <- act(model, s, action = a, fast = TRUE)
        s <- result$state_prime
        r <- result$reward
        rew <- rew + r * disc^(j - 1L)  
        
        
        if (trajectories) {
          trajectory[j, ] <-
            data.frame(
              episode = i,
              time = j - 1L,
              s = features2state(s_prev),
              a = a,
              r = r,
              s_prime = features2state(s)
            )
        }
        
        if (absorbing_states(model, s)) {
          if (trajectories) {
            trajectory <- trajectory[1:j, , drop = FALSE]
          }
          break
        }
      }
      
      list(
        action_cnt = action_cnt,
        state_cnt = state_cnt,
        reward = rew,
        trajectory = trajectory
      )
    }
    
    rew <- Reduce(c, lapply(sim, "[[", "reward"))
    rew <- unname(rew)
    
    the_trajectories <- NULL
    if (trajectories) {
      the_trajectories <- Reduce(rbind, lapply(sim, "[[", "trajectory"))
      the_trajectories$a <- normalize_action(the_trajectories$a, model)
    }
   
    # reduce the state counts by adding them all to a new fast map
    state_cnt <- fastmap()
    for (i in seq_along(sim)) {
      keys <- (sim[[i]]$state_cnt)$keys()
      joined_list <- unlist(state_cnt$mget(keys, missing = 0L)) + unlist((sim[[i]]$state_cnt)$mget(keys, missing = 0L))
      joined_list <- split(unname(joined_list), names(joined_list))
      state_cnt$mset(.list = joined_list)
    }
    state_cnt <- unlist(state_cnt$as_list(sort = TRUE))
     
     
    samp <- list(
      avg_reward = mean(rew, na.rm = TRUE),
      reward = rew,
      action_cnt = Reduce("+", lapply(sim, "[[", "action_cnt")),
      state_cnt = state_cnt, 
      trajectories = the_trajectories
    )
    
    samp$avg_episode_length = sum(samp$state_cnt) / n
    
    return(samp)
    
  }
