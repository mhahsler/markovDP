# Explanation for how to parametrize q(s,a):
# https://danieltakeshi.github.io/2016/10/31/going-deeper-into-reinforcement-learning-understanding-q-learning-and-linear-function-approximation/
# See file linear_function_approximation!


#' Solve MDPs with Temporal Differencing with Function Approximation
#'
#' Solve the MDP control problem using state-value approximation
#' by semi-gradient Sarsa (temporal differencing) for
#' episodic problems.
#'
#' 
#' ## Episodic Semi-gradient Sarsa
#'
#' The implementation follows the temporal difference algorithm 
#' episodic Semi-gradient Sarsa algorithm given in Sutton and Barto (2018).
#' 
#' ## Schedules
#' 
#' * epsilon schedule: `t` is increased by each processed episode.
#' * alpha schedule: `t` is increased by each processed episode.
#' 
#' @family solver
#' @family MDPTF
#'
#' @references
#' Sutton, Richard S., and Andrew G. Barto. 2018. Reinforcement Learning: An Introduction. Second. The MIT Press. [http://incompleteideas.net/book/the-book-2nd.html](http://incompleteideas.net/book/the-book-2nd.html).

#' @examples
#' # Example 1: A maze without walls. The step cost is 1. The start is top-left and
#' # the goal (+100 reward) is bottom-right.
#' # This is the ideal problem for a linear approximation of the Q-function
#' # using the x/y location as state features.
#'
#' m <- gw_maze_MDP(c(5, 5), start = "s(1,1)", goal = "s(5,5)")
#' 
#' # gridworlds have state labels om the format "s(row, col)" which can be
#' # automatically converted into state features used for approximation.
#' S(m)
#' get_state_features(m)
#'
#' # solve using linear state features (no transformation) 
#' set.seed(1000)
#' sol <- solve_MDP_APPROX(m, horizon = 1000, n = 100) 
#'
#' # approximation
#' sol$solution$q_approx_linear
#' 
#' gw_plot(sol)
#' gw_matrix(sol, what = "value")
#'
#' # the approximate value function can be visualized for states 
#' # with two features.
#' approx_V_plot(sol)
#'
#' # extracting approximate Q-values
#' approx_greedy_action(sol, "s(4,5)")
#' approx_Q_value(sol, "s(4,5)", "down")
#' approx_Q_value(sol)
#'
#' # extracting a greedy policy using the approximate Q-values
#' approx_greedy_policy(sol)
#'
#'
#' # Example 2: Stuart Russell's 3x4 Maze using linear basis approximation
#' # The wall and the -1 absorbing state make linear approximation
#' # using just the position directly more difficult.
#'
#' data(Maze)
#' gw_plot(Maze)
#'
#' # if no state features are specified, then they are constructed
#' # by parsing the state label of the form s(feature list).
#' set.seed(1000)
#' sol <- solve_MDP_APPROX(Maze, horizon = 100, n = 100,
#'                         alpha = schedule_exp(0.3, 0.01),
#'                         epsilon = schedule_exp(1, 0.1))
#' gw_plot(sol)
#' gw_matrix(sol, what = "value")
#' approx_V_plot(sol, res = 20)
#'
#'
#' # Example 3: Stuart Russell's 3x4 Maze using
#' #            order-1 Fourier basis for approximation and
#' #            1-step Sarsa
#'
#' set.seed(1000)
#' sol <- solve_MDP_APPROX(Maze, horizon = 100, n = 100,
#'                     alpha = schedule_exp(0.3, .01),
#'                     epsilon = schedule_exp(1, .1),
#'                     transformation = transformation_fourier_basis, 
#'                     order = 1
#'                     )
#'                     
#' gw_plot(sol)
#' gw_matrix(sol, what = "value")
#' approx_V_plot(sol, res = 20)
#'
#' # Example 4: Stuart Russell's 3x4 Maze using
#' #    order-1 Fourier basis for approximation
#' #    and eligibility traces: Sarsa(lambda)
#' 
#' set.seed(1000)
#'
#' ## TODO: The following example does not converge to 1!
#' data(Maze)
#' sol <- solve_MDP_APPROX(Maze, horizon = 100, n = 100,
#'                     alpha = schedule_exp(0.3, .01),
#'                     epsilon = schedule_exp(1, .1),
#'                     lambda = 0.1,
#'                     transformation = transformation_fourier_basis, 
#'                     order = 1
#'                     )
#'                     
#' gw_plot(sol)
#' gw_matrix(sol, what = "value")
#' approx_V_plot(sol, res = 20)
#'  
#' @inheritParams solve_MDP_TD
#' @param method string; one of the following solution methods: `'sarsa'`
#' @param lambda the trace-decay parameter for the an accumulating trace. If `lambda = 0`
#'      then 1-step Sarsa is used.
#' @param w an initial weight vector. By default a vector with 0s is used.
#' @param transformation a transformation function. See [transformation].
#' @param ... further parameters are passed on to the [transformation] function. 
#' @param state_features a matrix with one row per state with state features 
#'      to be used. If `NULL` then [get_state_features()] will be used to get 
#'      the state features stored in the model, or to construct 
#'      state features from state labels.
#'
#' @inherit solve_MDP return
#'
#' @export
solve_MDP_APPROX <-
  function(model,
           method = "sarsa",
           horizon = NULL,
           discount = NULL,
           alpha = schedule_exp(0.2, .1),
           epsilon = schedule_exp(1, .1),
           lambda = 0,
           n,
           state_features = NULL,
           transformation = transformation_linear_basis,
           w = NULL,
           ...,
           matrix = TRUE,
           continue = FALSE,
           progress = TRUE,
           verbose = FALSE) {
    
    if (lambda == 0)
      solve_MDP_APPROX_1_step(model,
                              method,
                              horizon,
                              discount,
                              alpha,
                              epsilon,
                              n,
                              state_features,
                              transformation,
                              w = NULL,
                              ...,
                              matrix = matrix,
                              continue = continue,
                              progress = progress,
                              verbose = verbose)
    else
      solve_MDP_APPROX_lambda(model,
                              method,
                              horizon,
                              discount,
                              alpha,
                              epsilon,
                              lambda,
                              n,
                              state_features,
                              transformation,
                              w = NULL,
                              ...,
                              matrix = matrix,
                              continue = continue,
                              progress = progress,
                              verbose = verbose)
      
  }
  
solve_MDP_APPROX_1_step <-
  function(model,
           method = "sarsa",
           horizon = NULL,
           discount = NULL,
           alpha = schedule_exp(0.2, .1),
           epsilon = schedule_exp(1, .1),
           n,
           state_features = NULL,
           transformation,
           w = NULL,
           ...,
           matrix = TRUE,
           continue = FALSE,
           progress = TRUE,
           verbose = FALSE) {
    
    if (!inherits(model, "MDPE"))
      stop("This model needs to be a MDP environment.")
    
    method <-
      match.arg(method, c("sarsa"))
    
    alpha_seq <- .schedule_to_sequence(alpha, n, continue, model)    
    epsilon_seq <- .schedule_to_sequence(epsilon, n, continue, model)
    
    # this code solve MDP and MDPTF
    model <- .prep_model(model, horizon, discount, matrix, verbose, progress)
    
    if (verbose)
      progress <- FALSE
    
    if (progress) {
      pb <- my_progress_bar(n + 1L, name = "solve_MDP")
      pb$tick(0)
    }
    
    if (!is.null(state_features))
      model$state_features <- state_features
    
    # Initialize approximation
    if (continue) {
      if (!is.null(w))
        stop("continue and w cannot be both specified!")
      q <- model$solution$q_approx_linear
      w <- w$w
      if (is.null(q))
        stop("Model does not contain an approximation function. Cannot continue!")
    } else {
      q <- q_approx_linear(model, transformation = transformation, ...)
      model$solution$q_approx_linear <- q
      w <- w %||% q$w 
    }
   
     
    # return unconverged result when interrupted
    on.exit({
      if (progress) {
        pb$tick(0)
        pb$terminate()
      }
      
      if (e < n)
        warning("Manual interupt: MDP solver stopped at episode ", e)
      
      if (verbose) {
        cat("\nTerminated at episode:", e, "\n")
      }
      
      q$w <- w
      
      model$solution <- list(
        method = method,
        alpha = alpha,
        epsilon = epsilon,
        n = n + ifelse(continue, model$solution$n, 0),
        converged = NA,
        q_approx_linear = q,
        policy = if (!is.null(S(model)))
          list(approx_greedy_policy(model, w = w))
        else
          NULL
      )
      return(model)
    })
    
    if (verbose) {
      cat("Running", method)
      cat("\nalpha:            ", show_schedule(alpha))
      cat("\nepsilon:          ", show_schedule(epsilon))
      cat("\nn:                ", n)
      if (continue) {
        cat("\ncont. with w:\n")
        print(w)
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
      alpha_val <- alpha_seq[e]
      
      # MDP: state is an id, MDPTF: state is a features
      s <- start(model)
      
      a <- approx_greedy_action(model, s, w, epsilon_val, as = "id")
      
      # loop steps in episode
      i <- 0L
      while (TRUE) {
        i <- i + 1L
        if ((i %% 100 == 0) && progress && !pb$finished)
          pb$tick(0)
        
        # take action a
        a_res <- act(model, s, a, fast = TRUE)
        s_prime <- a_res$state_prime
        r <- a_res$r
        
        # unavailable actions return a reward of -inf
        if (!is.finite(r))
          stop("Inf reward not supported!")
        
        s_features <- normalize_state_features(s, model)
        
        # reached terminal state?
        if (absorbing_states(model, state = s_prime)) {
          if (verbose > 1) {
            cat("\nTerminal state reached!")
            cat(
              sprintf(
                "\nStep %i - s=%s a=%s r=%.3f s'=%s - q(s,a)=%.3f error=%.3f\n",
                i,
                normalize_state_label(s, model),
                normalize_action_label(a, model),
                r,
                normalize_state_label(s_prime, model),
                q$f(s_features, a, w),
                r - q$f(s_features, a, w)
              )
            )
            print(w)
          }
          
          w <- w + alpha_val * (r - q$f(s_features, a, w)) * q$gradient(s_features, a, w)
          
          if (verbose > 1) {
            print(w)
            cat(sprintf("New q(s,a)=%.3f\n", q$f(s_features, a, w)))
          }
          
          break
        }
        
        s_prime_features <- normalize_state_features(s_prime, model)
        a_prime <- approx_greedy_action(model, s_prime, w, epsilon_val, as = "id")
        
        if (verbose > 1) {
          if (i == 1L) {
            cat("\n*** Episode", e, "***\n")
          }
          cat(
            sprintf(
              "\nStep %i - s=%s a=%s r=%.3f s'=%s a'=%s - q(s,a)=%.3f q(s',a')=%.3f error=%.3f\n",
              i,
              normalize_state_label(s, model),
              normalize_action_label(a, model),
              r,
              normalize_state_label(s_prime, model),
              normalize_action_label(a_prime, model),
              q$f(s_features, a, w),
              q$f(s_prime_features, a_prime, w),
              r + discount * q$f(s_prime_features, a_prime, w) - q$f(s_features, a, w)
            )
          )
          print(w)
        }
        
        # this is one-step Sarsa
        
        # note:  q$f(s_prime, a_prime, w) for absorbing states can be very off!
        td_error <- (r + discount * q$f(s_prime_features, a_prime, w) - q$f(s_features, a, w))
        if (td_error > 1e20)
          stop(
            "Temporal differencing error becomes too large which indicates instability. Reduce alpha."
          )
        
        w <- w + alpha_val * td_error * q$gradient(s_features, a, w)
        
        if (verbose > 1) {
          print(w)
          cat(sprintf("New q(s,a)=%.3f\n", q$f(s_features, a, w)))
        }
        
        s <- s_prime
        a <- a_prime
        
        if (i >= horizon)
          break
        
      }
    }
    
    # return via on.exit()
  }


solve_MDP_APPROX_lambda <-
  function(model,
           method = "sarsa",
           horizon = NULL,
           discount = NULL,
           alpha = schedule_exp(0.2, .1),
           epsilon = schedule_exp(1, .1),
           lambda = 0.1,
           n,
           state_features = NULL,
           transformation,
           w = NULL,
           ...,
           matrix = TRUE,
           continue = FALSE,
           progress = TRUE,
           verbose = FALSE) {
    
    if (!inherits(model, "MDPE"))
      stop("This model needs to be a MDP environment.")
    
    methods <- c("Sarsa", "GTD")
    method <- match.arg(tolower(method), tolower(methods))
    method_id <- match(tolower(method), tolower(methods))
    
    alpha_seq <- .schedule_to_sequence(alpha, n, continue, model)    
    epsilon_seq <- .schedule_to_sequence(epsilon, n, continue, model)
    
    # this code solve MDP and MDPTF
    model <- .prep_model(model, horizon, discount, matrix, verbose, progress)
    
    if (verbose)
      progress <- FALSE
    
    if (progress) {
      pb <- my_progress_bar(n + 1L, name = "solve_MDP")
      pb$tick(0)
    }
    
    if (!is.null(state_features))
      model$state_features <- state_features
    
    # Initialize approximation
    if (continue) {
      if (!is.null(w))
        stop("continue and w cannot be both specified!")
      qf <- model$solution$q_approx_linear
      w <- qf$w
      if (is.null(q))
        stop("Model does not contain an approximation function. Cannot continue!")
    } else {
      qf <- q_approx_linear(model, transformation = transformation, ...)
      model$solution$q_approx_linear <- qf
      w <- w %||% qf$w 
    }
    
    
    # return unconverged result when interrupted
    on.exit({
      if (progress) {
        pb$tick(0)
        pb$terminate()
      }
      
      if (e < n)
        warning("Manual interupt: MDP solver stopped at episode ", e)
      
      if (verbose) {
        cat("\nTerminated at episode:", e, "\n")
      }
      
      qf$w <- w
      
      model$solution <- list(
        method = method,
        alpha = alpha,
        epsilon = epsilon,
        lambda = lambda,
        n = n + ifelse(continue, model$solution$n, 0),  
        q_approx_linear = qf,
        converged = NA,
        policy = if (!is.null(S(model)))
          list(approx_greedy_policy(model, w = w))
        else
          NULL
      )
      return(model)
    })
    
    if (verbose) {
      cat("Running", method)
      cat("\nalpha:            ", show_schedule(alpha))
      cat("\nepsilon:          ", show_schedule(epsilon))
      cat("\nlambda:           ", lambda)
      cat("\nn:                ", n)
      if (continue) {
        cat("\ncont. with w:\n")
        print(w)
      }
    }
    
    A <- A(model)
    gamma <- model$discount
    
    # loop for each episodes
    e <- 0L
    while (e < n) {
      e <- e + 1L
      if (progress)
        pb$tick()
      
      epsilon_val <- epsilon_seq[e]
      alpha_val <- alpha_seq[e]
      
      # Initialize s and choose first action
      # MDP: state is an id, MDPTF: state is a features
      s <- start(model)
      a <- approx_greedy_action(model, s, w, epsilon_val, as = "id")
      
      x <- qf$x(normalize_state_features(s, model), a)
      
      # initialize z to all 0s
      z <- qf$w
      z[] <- 0
      
      q_old <- 0
      
      # loop steps in episode
      i <- 0L
      while (TRUE) {
        i <- i + 1L
        if ((i %% 100 == 0) && progress && !pb$finished)
          pb$tick(0)
        
        # take action a and observe r and s'
        a_res <- act(model, s, a, fast = TRUE)
        s_prime <- a_res$state_prime
        r <- a_res$r
        
        # choose a'
        a_prime <- approx_greedy_action(model, s_prime, w, epsilon_val, as = "id")
        
        # get x'
        x_prime <- qf$x(normalize_state_features(s_prime, model), a_prime)
        
        q <- sum(x * w)
        q_prime <- sum(x_prime * w)
        delta <- r + gamma * q_prime - q
        
        if (delta > 1e20)
          stop(
            "Temporal differencing error becomes too large which indicates instability. Reduce alpha or lambda."
          )
        
        z <- gamma * lambda * z + (1 - alpha_val * gamma * lambda * sum(z*x)) * x
        
        # Sarsa(lambda) .. on-policy true online
        if (method_id == 1L) {
          w <- w + alpha_val * (delta + q - q_old) * z - alpha_val * (q - q_old) * x
        }
        
        # GTD(lambda)
        else if (method_id == 2) {
          w <- w + alpha_val * delta * z - alpha_val * gamma * (1 - lambda) * sum(z*q) * x
        }
        
        # GQ(lambda)
        # average is only one since the target policy is not soft.
        #x_bar <- qf$x(normalize_state_features(s, model), approx_greedy_action(model, s, w, epsilon = 0))
        #delta_a <- r + gamma * sum(w*x_bar - sum(w*x))
        #w <- w + alpha * delta_a * z - alpha * gamma * (1 - lambda) * (sum(z*q)) * x_bar
        #
        
        if (verbose > 1) {
          if (i == 1L) {
            cat("\n*** Episode", e, "***\n")
          }
          cat(
            sprintf(
              "\nStep %i - s=%s a=%s r=%.3f s'=%s a'=%s - q(s,a)=%.3f q(s',a')=%.3f delta=%.3f\n",
              i,
              normalize_state_label(s, model),
              normalize_action_label(a, model),
              r,
              normalize_state_label(s_prime, model),
              normalize_action_label(a_prime, model),
              q,
              q_prime,
              delta
            )
          )
          cat("w:\n")
          print(w)
          cat("z:\n")
          print(z)
        }
        
        if (i >= horizon)
          break
        
        if(absorbing_states(model, state = s_prime))
          break
        
        q_old <- q_prime
        x <- x_prime
        s <- s_prime
        a <- a_prime
        
      }
    }
    
    # return via on.exit()
  }


#' @rdname solve_MDP_APPROX
#' @param state a state (index or name)
#' @param action an action (index or name)
#' @param w a weight vector
#' @param as character; specifies the desired output format (see [normalize_action()])
#' @export
approx_Q_value <- function(model,
                           state = NULL,
                           action = NULL,
                           w = NULL) {
  action <- action %||% A(model)
  
  state <- state %||% S(model)
  if (is.null(state))
    stop("a state needs to be specified!")
  
  w <- w %||% model$solution$q_approx_linear$w
  if (is.null(w))
    stop("weight vector w is missing.")
  
  if(!is.matrix(state))
    state <- normalize_state_features(state, model)
  
  if (length(action) > 1L)
    return(sapply(
      action,
      FUN = function(a)
        approx_Q_value(model, state, a, w)
    ))
  
  if (nrow(state) == 1L)
    model$solution$q_approx_linear$f(drop(state), action, w)
  else
    apply(state, MARGIN = 1, model$solution$q_approx_linear$f, action, w)
}

#' @rdname solve_MDP_APPROX
#' @export
approx_greedy_action <- function(model,
                                 state,
                                 w = NULL,
                                 epsilon = 0,
                                 as = "factor") {
  if (epsilon == 0 || runif(1) > epsilon) {
    a <- which.max.random(approx_Q_value(model, state, w = w))
  } else {
    a <- sample.int(length(A(model)), 1L)
  }
  
  normalize_action(a, model, as)
}

#' @rdname solve_MDP_APPROX
#' @export
approx_greedy_policy <- function(model, w = NULL) {
  w <- w %||% model$solution$w
  S <- S(model)
  
  if (is.null(S))
    stop("Policy is only available for a specified finite state space!")
  
  qs <- approx_Q_value(model, w = w)
  
  data.frame(
    state = S,
    V = apply(qs, MARGIN = 1, max),
    action = normalize_action(apply(qs, MARGIN = 1, which.max.random), model),
    row.names = NULL
  )
}


# internal: transpose and reorder rows for a proper image
#' @importFrom graphics axis contour image
pimage <- function (x1,
                    x2,
                    z,
                    col = hcl.colors(12, "YlOrRd", rev = TRUE),
                    image = TRUE,
                    contour = TRUE,
                    axes = TRUE,
                    ...) {
  
  z <- t(z)[, rev(seq_len(nrow(z))), drop = FALSE]
  
  if (image) {
    image(x2, x1, z, axes = FALSE, col = col, ...)
    
    if (contour)
      contour(x2, x1, z, axes = FALSE, add = TRUE)
  } else
    contour(x2, x1, z, axes = FALSE, ...)
  
  box()
  
  if (axes) {
    p <- pretty(x2)
    axis(1L,
         at = seq(min(p), max(p), along.with = p),
         labels = p)
    p <- pretty(x1)
    axis(2L,
         at = seq(min(p), max(p), along.with = p),
         labels = p)
  }
}

#' @rdname solve_MDP_APPROX
#' @param image,contour logical; include the false color image or the 
#'        contours in the plot?
#' @param main title for the plot.
#' @param res resolution as the number of values sampled from each feature.
#' @param col colors for the passed on to [`image()`].
#' @param min,max numeric vectors with minimum/maximum values for each feature
#'        in the state feature representation.
#' @export
approx_V_plot <- function(model,
                          min = NULL,
                          max = NULL,
                          w = NULL,
                          res = 25,
                          col = hcl.colors(res, "YlOrRd", rev = TRUE),
                          image = TRUE,
                          contour = TRUE,
                          main = NULL,
                          ...) {
  if (is.null(main)) {
    main <- paste("Approx. Value Function:",
                  model$name,
                  paste0("(", model$solution$method, ")"))
  }
  
  rng <- get_state_feature_range(model, min, max)
  min <- rng[1, ]
  max <- rng[2, ]
  dim_s <- length(max)
  
  if (dim_s != 2)
    stop("This visualization only works for states with 2 features!")
  
  x1 <- seq(min[1], max[1], length.out = res)
  x2 <- seq(min[2], max[2], length.out = res)
  
  V_val <- Vectorize(function(x1, x2)
    max(approx_Q_value(model, state = s(x1, x2), w = w)))
  V <- outer(x1, x2, V_val)
  
  pimage(x1, x2, V, col, image, contour, main = main, ...)
}