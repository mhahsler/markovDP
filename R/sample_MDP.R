#' Sample Trajectories from an MDP
#'
#' Sample trajectories through an MDP. The start state for each
#' trajectory is randomly chosen using the specified belief. The belief is used to choose actions
#' from an epsilon-greedy policy and then update the state.
#'
#' The default is a
#' faster C++ implementation (`engine = 'cpp'`).
#' A native R implementation is available (`engine = 'r'`).
#'
#' Both implementations support parallel execution using the package
#' \pkg{foreach}. To enable parallel execution, a parallel backend like
#' \pkg{doparallel} needs to be available needs to be registered (see
#' [doParallel::registerDoParallel()]).
#' Note that small samples are slower using parallelization. Therefore, C++ simulations
#' with n * horizon less than 100,000 are always executed using a single worker.
#' @family MDP
#' @importFrom stats runif
#'
#' @param model an MDP model.
#' @param n number of trajectories.
#' @param start probability distribution over the states for choosing the
#'  starting states for the trajectories. Defaults to "uniform".
#' @param horizon epochs end once an absorbing state is reached or after
#'  the maximal number of epochs specified via `horizon`. If `NULL` then the
#'  horizon for the model is used.
#' @param epsilon the probability of random actions  for using an epsilon-greedy policy.
#'  Default for solved models is 0 and for unsolved model 1.
#' @param engine `'cpp'` or `'r'` to perform simulation using a faster C++
#'  or a native R implementation `NULL` uses the C++ implementation unless the transition model or
#'  the reward are specified as R functions (which are slow in C++).
#' @param trajectories logical; return the complete trajectories.
#' @param delta_horizon precision used to determine the horizon for infinite-horizon problems.
#' @param exploring_starts logical; randomly sample a start/action combination to
#'  start the episode from.
#' @param progress show a progress bar?
#' @param verbose report used parameters
#' @param ... further arguments are ignored.
#' @return A list with elements:
#'  * `avg_reward`: The average discounted reward.
#'  * `reward`: Reward for each trajectory.
#'  * `action_cnt`: Action counts.
#'  * `state_cnt`: State counts.
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
#' data(Maze)
#'
#' # solve the MDP for 5 epochs and no discounting
#' sol <- solve_MDP(Maze, discount = 1)
#' sol
#'
#' # V in the policy is and estimate of the state values when following the optimal policy.
#' policy(sol)
#' gw_matrix(sol, what = "action")
#'
#' ## Example 1: simulate 100 trajectories following the policy,
#' #             only the final belief state is returned
#' sim <- sample_MDP(sol, n = 100, horizon = 10, verbose = TRUE)
#' sim
#'
#' # Note that all simulations for this model start at s_1 and that the simulated avg. reward
#' # is therefore an estimate to the value function for the start state s_1.
#' policy(sol)[1, ]
#'
#' # Calculate proportion of actions taken in the simulation
#' round_stochastic(sim$action_cnt / sum(sim$action_cnt), 2)
#'
#' # reward distribution
#' hist(sim$reward)
#'
#' ## Example 2: simulate starting following a uniform distribution over all
#' #             states and return all trajectories
#' sim <- sample_MDP(sol,
#'   n = 100, start = "uniform", horizon = 10,
#'   trajectories = TRUE
#' )
#' head(sim$trajectories)
#'
#' # how often was each state visited?
#' table(sim$trajectories$s)
#' @export
sample_MDP <-
  function(model,
           n = 100,
           start = NULL,
           horizon = NULL,
           epsilon = NULL,
           exploring_starts = FALSE,
           delta_horizon = 1e-3,
           trajectories = FALSE,
           engine = NULL,
           progress = TRUE,
           verbose = FALSE,
           ...) {
    .nodots(...)
    
    if (is.null(engine)) {
      if (is.function(model$transition_prob) ||
          is.function(model$reward))
        engine <- "r"
      else 
        engine <- "cpp"
    }
    
    engine <- match.arg(tolower(engine), c("cpp", "r"))
    
    solved <- is_solved_MDP(model)
    n <- as.integer(n)
    
    # exploring starts: uniform start + first action is random
    if (exploring_starts) {
      if (!is.null(start))
        warning("start cannot be specified for exploring starts. Using 'uniform'!")
      start <- "uniform"
    }
    
    start <- start_vector(model, start = start, sparse = FALSE)
    
    if (is.null(horizon)) {
      horizon <- model$horizon
    }
    if (is.null(horizon) || is.infinite(horizon)) {
      if (!is.null(model$discount) && model$discount < 1) {
        # find a horizon that approximates the reward using
        # discount^horizon * max_abs_R <= 0.001
        max_abs_R <- max(abs(.reward_range(model)))
        horizon <-
          ceiling(log(delta_horizon / max_abs_R) / log(model$discount))
      } else { 
        warning("Simulation for undiscounted problems need a finite simulation horizon.\n",
                "Using 10000.")
        horizon <- 10000
      }
    }
    horizon <- as.integer(horizon)
    
    if (is.null(epsilon)) {
      if (!solved) {
        epsilon <- 1
      } else {
        epsilon <- 0
      }
    }
    
    if (!solved && epsilon != 1) {
      stop("epsilon has to be 1 for unsolved models.")
    }
    
    disc <- model$discount
    if (is.null(disc)) {
      disc <- 1
    }
    
    if (engine == "cpp") {
      if (foreach::getDoParWorkers() == 1 || n * horizon < 100000) {
        return(
          sample_MDP_cpp(
            model,
            n,
            start,
            horizon,
            disc,
            trajectories,
            epsilon,
            exploring_starts,
            verbose = verbose
          )
        )
      }
      
      ns <- foreach_split(n)
      
      if (verbose) {
        cat("Sampling MDP trajectories.\n")
        cat("- engine: cpp \n")
        cat("- horizon:", horizon, "\n")
        cat("- n:", n, "- parallel workers:", length(ns), "\n")
        cat("- epsilon:", epsilon, "\n")
        cat("- discount factor:", disc, "\n")
        cat("\n")
      }
      
      w <-
        NULL # to shut up the warning for the foreach counter variable
      
      sim <- foreach(w = 1:length(ns)) %dopar%
        sample_MDP_cpp(
          model,
          ns[w],
          start,
          horizon,
          disc,
          trajectories,
          epsilon,
          exploring_starts,
          verbose = FALSE
        )
      
      # adjust the episode number for parallel processing
      episode_add <- cumsum(c(0L, ns))
      for (i in seq_along(sim)) {
        sim[[i]]$trajectories$episode <-
          sim[[i]]$trajectories$episode + episode_add[i]
      }
      
      rew <- Reduce(c, lapply(sim, "[[", "reward"))
      
      return(
        list(
          avg_reward = mean(rew, na.rm = TRUE),
          reward = rew,
          action_cnt = Reduce("+", lapply(sim, "[[", "action_cnt")),
          state_cnt = Reduce("+", lapply(sim, "[[", "state_cnt")),
          trajectories = Reduce(rbind, lapply(sim, "[[", "trajectories"))
        )
      )
    }
    
    # R implementation starts here ##############
    
    states <- as.character(S(model))
    n_states <- length(states)
    states_absorbing <- absorbing_states(model, sparse = "index")
    actions <- as.character(A(model))
    
    # for easier access
    pol <-
      lapply(
        model$solution$policy,
        FUN = function(p) {
          structure(p$action, names = p$state)
        }
      )
    
    if (verbose) {
      cat("Sampling MDP trajectories.\n")
      cat("- engine:", engine, "\n")
      cat("- horizon:", horizon, "\n")
      cat("- exploring starts:", exploring_starts, "\n")
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
      pb <- my_progress_bar(n, name = "sample_MDP")
    
    #warning("Debug mode on!!!")
    #sim <- for(i in 1:n){
    sim <- foreach(i = 1:n) %dopar% {
      if (progress)
        pb$tick()
      
      # find a initial state
      s <- sample.int(length(states), 1L, prob = start)
      
      action_cnt <- rep(0L, length(actions))
      names(action_cnt) <- actions
      state_cnt <- rep(0L, length(states))
      names(state_cnt) <- states
      rew <- 0
      
      if (trajectories) {
        trajectory <- data.frame(
          episode = rep(NA_integer_, horizon),
          time = rep(NA_integer_, horizon),
          s = NA_integer_,
          a = NA_integer_,
          r = NA_real_,
          s_prime = NA_integer_
        )
      } else {
        trajectory <- NULL
      }
      
      for (j in seq_len(horizon)) {
        if (exploring_starts && j == 1L) {
          # choose first action randomly
          a <- sample.int(length(actions), 1L, replace = TRUE)
        } else {
          # epsilon soft policy
          if (runif(1) < epsilon) {
            a <- sample.int(length(actions), 1L, replace = TRUE)
          } else {
            a <- pol[[.get_pol_index(model, j)]][s]
          }
        }
        
        action_cnt[a] <- action_cnt[a] + 1L
        state_cnt[s] <- state_cnt[s] + 1L
        
        s_prev <- s
        s <-
          sample.int(length(states), 1L, prob = transition_matrix(model, a, s, sparse = FALSE))
        
        # rew <- rew + rew_m[[a]][[s_prev]][s] * disc ^ (j - 1L)
        # MDPs have no observation!
        r <- reward_matrix(model, a, s_prev, s)
        rew <- rew + r * disc ^ (j - 1L)
        
        if (trajectories) {
          trajectory[j, ] <-
            data.frame(
              episode = i,
              time = j - 1L,
              s = s_prev,
              a = a,
              r = r,
              s_prime = s
            )
        }
        
        if (s %in% states_absorbing) {
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
      the_trajectories$s <- .normalize_state(the_trajectories$s, model)
      the_trajectories$a <- .normalize_action(the_trajectories$a, model)
      the_trajectories$s_prime <- .normalize_state(the_trajectories$s_prime, model)
    }
    
    list(
      avg_reward = mean(rew, na.rm = TRUE),
      reward = rew,
      action_cnt = Reduce("+", lapply(sim, "[[", "action_cnt")),
      state_cnt = Reduce("+", lapply(sim, "[[", "state_cnt")),
      trajectories = the_trajectories
    )
  }
