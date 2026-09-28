#' Solve MDPs with Policy Gradient Methods
#'
#' Not implemented yet: Solve the MDP control problem using
#' a parameterized policy and policy gradient methods.
#' The implemented method one-step actor-critic control and
#' actor-critic control with eligibility traces.
#'
#'
#' @family solver
#' @family MDPSample
#'
#' @references
#' Sutton, Richard S., and Andrew G. Barto. 2018. Reinforcement Learning: An Introduction. Second. The MIT Press. [http://incompleteideas.net/book/the-book-2nd.html](http://incompleteideas.net/book/the-book-2nd.html).
#'
#' @examples
#' # EXAMPLE
#' @inheritParams solve_MDP_TD
#' @param method string; one of the following solution methods: `'sarsa'`
#' @param lambda the trace-decay parameter for the an accumulating trace. If `lambda = 0`
#'      then 1-step Sarsa is used.
#' @param w an initial weight vector. By default a vector with 0s is used.
#' @param theta parameter...
#' @param alpha_actor,alpha_critic alpha schedules
#'
#' @inherit solve_MDP return
#'
#' @export
solve_MDP_PG <-
  function(model,
           method = "actor-critic",
           horizon = NULL,
           discount = NULL,
           alpha_actor = schedule_exp(0.2, .1),
           alpha_critic = schedule_exp(0.2, .1),
           epsilon = schedule_exp(1, .1),
           lambda = 0,
           n,
           w = NULL,
           theta = NULL,
           ...,
           matrix = TRUE,
           continue = FALSE,
           progress = TRUE,
           verbose = FALSE) {
    .nodots(...)
    
    if (!inherits(model, "MDP"))
      stop("This model needs to be an MDP.")
    
    method <-
      match.arg(method, c("actor-critic"))
    
    warning("remove epsilon?")
    
    if (lambda == 0)
      solve_MDP_PG_AC_1_step(
        model = model,
        method = method,
        horizon = horizon,
        discount = discount,
        alpha_actor = alpha_actor,
        alpha_critic = alpha_critic,
        epsilon = epsilon,
        n = n,
        w_value = w,
        w_policy = theta,
        ...,
        matrix = matrix,
        continue = continue,
        progress = progress,
        verbose = verbose
      )
    else
      solve_MDP_PG_AC_lambda(
        model = model,
        method = method,
        horizon = horizon,
        discount = discount,
        alpha_actor = alpha_actor,
        alpha_critic = alpha_critic,
        epsilon = epsilon,
        lambda = lambda,
        n = n,
        w = w,
        theta = theta,
        ...,
        matrix = matrix,
        continue = continue,
        progress = progress,
        verbose = verbose
      )
    
  }

solve_MDP_PG_AC_1_step <-
  function(model,
           method = "sarsa",
           horizon = NULL,
           discount = NULL,
           alpha_actor = schedule_exp(0.2, .1),
           alpha_critic = schedule_exp(0.2, .1),
           epsilon = schedule_exp(1, .1),
           n,
           state_features = NULL,
           transformation_value_approx = transformation_linear_basis,
           w_value = NULL,
           transformation_policy_approx = transformation_linear_basis,
           w_policy = NULL,
           ...,
           matrix = TRUE,
           continue = FALSE,
           progress = TRUE,
           verbose = FALSE) {
    .nodots(...)
    stop("DO IT!")
    
    alpha_actor_seq <- .schedule_to_sequence(alpha_actor, n, continue, model)
    alpha_critic_seq <- .schedule_to_sequence(alpha_critic, n, continue, model)
    epsilon_seq <- .schedule_to_sequence(epsilon, n, continue, model)
    
    # this code solve MDP and MDPSample
    model <- .prep_model(model, horizon, discount, matrix, verbose, progress)
    
    if (verbose)
      progress <- FALSE
    
    if (progress) {
      pb <- my_progress_bar(n + 1L, name = "solve_MDP")
      pb$tick(0)
    }
    
    if (!is.null(state_features))
      model$state_features <- state_features
    
    
    # Initialize approximations
    if (continue) {
      if (!is.null(w_value) || !is.null(w_policy))
        stop("continue and w values cannot be both specified!")
      
      value_approx <- model$solution$value_approx
      policy_approx <- model$solution$policy_approx
      if (is.null(value_approx))
        stop("Model does not contain a value approximation function. Cannot continue!")
      if (is.null(policy_approx))
        stop("Model does not contain a policy approximation function. Cannot continue!")
      
      w_value <- value_approx$w
      w_policy <- policy_approx$w
      
    } else {
      value_approx <- v_approx_linear(model, transformation = transformation_value_approx, ...)
      policy_approx <- pi_approx_linear(model, transformation = transformation_policy_approx, ...)
      model$solution$value_approx <- value_approx
      model$solution$policy_approx <- policy_approx
      w_value <- w_value %||% value_approx$w
      w_policy <- w_policy %||% policy_approx$w
    }
    
    # return unconverged result when interrupted
    on.exit({
      if (progress) {
        pb$tick(0)
        pb$terminate()
      }
      
      if (e < n)
        warning("Manual interupt: MDP solver stopped at episode ", e)
      
      if (verbose)
        cat("\nTerminated at episode:", e, "\n")
      
      value_approx$w <- w_value
      policy_approx$w <- w_policy
      
      model$solution <- list(
        method = method,
        alpha_actor = alpha_actor,
        alpha_critic = alpha_critic,
        epsilon = epsilon,
        n = n + ifelse(continue, model$solution$n, 0),
        converged = NA,
        value_approx = value_approx,
        policy_approx = policy_approx,
        policy = if (!is.null(S(model)))
          list(approx_greedy_policy(model, w = w_policy))
        else
          NULL
      )
      return(model)
    })
    
    if (verbose) {
      cat("Running", method)
      cat("\nalpha_actor:      ", show_schedule(alpha_actor))
      cat("\nalpha_critic:     ", show_schedule(alpha_critic))
      cat("\nepsilon:          ", show_schedule(epsilon))
      cat("\nn:                ", n)
      if (continue) {
        cat("\ncont. with w_value:\n")
        print(w_value)
        cat("\ncont. with w_policy:\n")
        print(w_policy)
      }
    }
    
    A <- A(model)
    discount <- model$discount
    
    
    # loop episodes
    e <- 0L
    while (e < n) {
      e <- e + 1L
      if (progress)
        pb$tick()
      
      epsilon_val <- epsilon_seq[e]
      alpha_actor_val <- alpha_actor_seq[e]
      alpha_critic_val <- alpha_critic_seq[e]
      
      # s is an ID for MDPModel, and state features for MDPSample!
      s <- start(model)
      I <- 1
     
      # loop steps in episode
      i <- 0L
      while (TRUE) {
        i <- i + 1L
        if ((i %% 100 == 0) && progress && !pb$finished)
          pb$tick(0)
        
        # greedy action
        a <- which.max(approx_value(model, policy_approx, w_policy))
        
        # take action a and observe r and s'
        a_res <- act(model, s, a, fast = TRUE)
        s_prime <- a_res$state_prime
        r <- a_res$reward
        
        # calculate the error delta
        delta <- r + 
          discount*approx_value(value_approx, s_prime, w = w_value, model = model) -
          approx_value(value_approx, s, w = w_value, model = model)
        
        # update w_value
        
        # update w_policy
        
        
        # end episode?
        if (i >= horizon)
          break
        
        if(absorbing_states(model, state = s_prime))
          break
        
        I <- I * discount
        s <- s_prime
      }
    } 
    
    #     # return via on.exit()
  }

    
    solve_MDP_PG_AC_lambda <-
      function(model,
               method = "sarsa",
               horizon = NULL,
               discount = NULL,
               alpha_actor = schedule_exp(0.2, .1),
               alpha_critic = schedule_exp(0.2, .1),
               epsilon = schedule_exp(1, .1),
               lambda = 0.1,
               n,
               w = NULL,
               theta = NULL,
               ...,
               matrix = TRUE,
               continue = FALSE,
               progress = TRUE,
               verbose = FALSE) {
        .nodots(...)
        
        stop("DO IT!")
        #
        #     if (!inherits(model, "MDP"))
        #       stop("This model needs to be an MDP environment.")
        #
        #     methods <- c("Sarsa", "GTD")
        #     method <- match.arg(tolower(method), tolower(methods))
        #     method_id <- match(tolower(method), tolower(methods))
        #
        #     ### alpha/epsilon func
        #     alpha_seq <- NULL
        #     if (is.function(alpha))
        #       alpha_seq <- alpha(seq(1, n))
        #     else if (length(alpha) == n)
        #       alpha_seq <- alpha
        #
        #     epsilon_seq <- NULL
        #     if (is.function(epsilon))
        #       epsilon_seq <- epsilon(seq(1, n))
        #     else if (length(epsilon) == n)
        #       epsilon_seq <- epsilon
        #     ###
        #
        #
        #     qf <- model$approx_Q_function
        #     if (is.null(qf$w_init) || is.null(qf$f) || is.null(qf$gradient))
        #       stop(
        #         "Approx q-function f, the gradient, or w_init are missing in the q_function element in the model!"
        #       )
        #
        #     # this code solve MDP and MDPSample
        #     model <- .prep_model(model, horizon, discount, matrix, verbose, progress)
        #
        #     if (verbose)
        #       progress <- FALSE
        #
        #     if (progress) {
        #       pb <- my_progress_bar(n + 1L, name = "solve_MDP")
        #       pb$tick(0)
        #     }
        #
        #     # Initialize w
        #     if (continue) {
        #       if (!is.null(w))
        #         stop("continue and w cannot be both specified!")
        #       w <- model$solution$w
        #     } else {
        #       w <- w %||% qf$w_init
        #     }
        #
        #     # return unconverged result when interrupted
        #     on.exit({
        #       if (progress) {
        #         pb$tick(0)
        #         pb$terminate()
        #       }
        #
        #       if (e < n)
        #         warning("Manual interupt: MDP solver stopped at episode ", e)
        #
        #       if (verbose) {
        #         cat("\nTerminated at episode:", e, "\n")
        #       }
        #
        #       model$solution <- list(
        #         method = method,
        #         alpha = alpha,
        #         epsilon = epsilon,
        #         lambda = lambda,
        #         n = n,
        #         w = w,
        #         converged = NA,
        #         policy = if (!is.null(S(model)))
        #           list(approx_greedy_policy(model, w = w))
        #         else
        #           NULL
        #       )
        #       return(model)
        #     })
        #
        #     if (verbose) {
        #       cat("Running", method)
        #       cat("\nalpha:            ", show_schedule(alpha))
        #       cat("\nepsilon:          ", show_schedule(epsilon))
        #       cat("\nlambda:           ", lambda)
        #       cat("\nn:                ", n)
        #       if (continue) {
        #         cat("\ncont. with w:\n")
        #         print(w)
        #       }
        #     }
        #
        #     A <- A(model)
        #     gamma <- model$discount
        #
        #     # loop for each episodes
        #     e <- 0L
        #     while (e < n) {
        #       e <- e + 1L
        #       if (progress)
        #         pb$tick()
        #
        #
        #       ### alpha/epsilon func
        #       if (!is.null(epsilon_seq))
        #         epsilon <- epsilon_seq[e]
        #       if (!is.null(alpha_seq))
        #         alpha <- alpha_seq[e]
        #       ###
        #
        #       # Initialize s and choose first action
        #       # MDP: state is an id, MDPSample: state is a features
        #       s <- start(model)
        #       a <- approx_greedy_action(model, s, w, epsilon)
        #
        #       x <- qf$x(normalize_state_features(s, model), a)
        #       z <- qf$w_init
        #       q_old <- 0
        #
        #       # loop steps in episode
        #       i <- 0L
        #       while (TRUE) {
        #         i <- i + 1L
        #         if ((i %% 100 == 0) && progress && !pb$finished)
        #           pb$tick(0)
        #
        #         # take action a and observe r and s'
        #         a_res <- act(model, s, a, fast = TRUE)
        #         s_prime <- a_res$state_prime
        #         r <- a_res$r
        #
        #         # choose a'
        #         a_prime <- approx_greedy_action(model, s_prime, w, epsilon)
        #
        #         # get x'
        #         x_prime <- qf$x(normalize_state_features(s_prime, model), a_prime)
        #
        #         q <- sum(x * w)
        #         q_prime <- sum(x_prime * w)
        #         delta <- r + gamma * q_prime - q
        #
        #         if (delta > 1e20)
        #           stop(
        #             "Temporal differencing error becomes too large which indicates instability. Reduce alpha or lambda."
        #           )
        #
        #         z <- gamma * lambda * z + (1 - alpha + gamma * lambda * sum(z *
        #                                                                       x)) * x
        #
        #         # Sarsa(lambda) .. on-policy true online
        #         if (method_id == 1L) {
        #           w <- w + alpha * (delta + q - q_old) * z - alpha * (q - q_old) * x
        #         }
        #
        #         # GTD(lambda)
        #         else if (method_id == 2) {
        #           w <- w + alpha * delta * z - alpha * gamma * (1 - lambda) * sum(z * q) * x
        #         }
        #
        #         # GQ(lambda)
        #         # average is only one since the target policy is not soft.
        #         #x_bar <- qf$x(normalize_state_features(s, model), approx_greedy_action(model, s, w, epsilon = 0))
        #         #delta_a <- r + gamma * sum(w*x_bar - sum(w*x))
        #         #w <- w + alpha * delta_a * z - alpha * gamma * (1 - lambda) * (sum(z*q)) * x_bar
        #         #
        #
        #         if (verbose > 1) {
        #           if (i == 1L) {
        #             cat("\n*** Episode", e, "***\n")
        #           }
        #           cat(
        #             sprintf(
        #               "\nStep %i - s=%s a=%s r=%.3f s'=%s a'=%s - q(s,a)=%.3f q(s',a')=%.3f delta=%.3f\n",
        #               i,
        #               normalize_state_label(s, model),
        #               normalize_action_label(a, model),
        #               r,
        #               normalize_state_label(s_prime, model),
        #               normalize_action_label(a_prime, model),
        #               q,
        #               q_prime,
        #               delta
        #             )
        #           )
        #           # print(w)
        #           print(z)
        #         }
        #
        #         if (i >= horizon)
        #           break
        #
        #         if (absorbing_states(model, state = s_prime))
        #           break
        #
        #         q_old <- q_prime
        #         x <- x_prime
        #         s <- s_prime
        #         a <- a_prime
        #
        #       }
        #     }
        #
        #     # return via on.exit()
      }
    
    solve_MDP_PG_AC_lambda_infinite_horizon <-
      function(model,
               method = "sarsa",
               horizon = NULL,
               discount = NULL,
               alpha_actor = schedule_exp(0.2, .1),
               alpha_critic = schedule_exp(0.2, .1),
               epsilon = schedule_exp(1, .1),
               lambda = 0.1,
               n,
               w = NULL,
               theta = NULL,
               ...,
               matrix = TRUE,
               continue = FALSE,
               progress = TRUE,
               verbose = FALSE) {
        .nodots(...)
        
        stop("DO IT!")
      }
    
