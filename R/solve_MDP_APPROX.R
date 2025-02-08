# Explanation for how to parametrize q(s,a):
# https://danieltakeshi.github.io/2016/10/31/going-deeper-into-reinforcement-learning-understanding-q-learning-and-linear-function-approximation/

#' Episodic Semi-gradient Sarsa with Linear Function Approximation
#'
#' MDP control using state-value approximation. Semi-gradient Sarsa for
#' episodic problems.
#'
#' ## Linear Approximation
#' The state-action value function is approximated by
#' \deqn{\hat{q}(s,a) = \boldsymbol{w}^\top\phi(s,a),}
#'
#' where \eqn{\boldsymbol{w} \in \mathbb{R}^n} is a weight vector
#' and \eqn{\phi: S \times A \rightarrow \mathbb{R}^n}  is a
#' feature function that
#' maps each state-action pair to a feature vector.
#' The gradient of the state-action function is
#' \deqn{\nabla \hat{q}(s,a,\boldsymbol{w}) = \phi(s,a).}
#'
#' ## State-action Feature Vector Construction
#'
#' For a small number of actions, we can
#' follow the construction described by Geramifard et al (2013)
#' which uses a state feature function \eqn{\phi: S \rightarrow \mathbb{R}^{m}}
#' to construct the complete state-action feature vector.
#' Here, we also add an intercept term.
#' The state-action feature vector has length \eqn{1 + |A| \times m}.
#' It has the intercept and then one component for each action. All these components
#' are set to zero and only the active action component is set to \eqn{\phi(s)},
#' where \eqn{s} is the current state.
#' For example, for the state feature
#' vector \eqn{\phi(s) = (3,4)} and action \eqn{a=2} out of three possible
#' actions \eqn{A = \{1, 2, 3\}}, the complete
#' state-action feature vector is \eqn{\phi(s,a) = (1,0,0,3,4,0,0)}.
#' The leading 1 is for the intercept and the zeros represent the two not
#' chosen actions.
#'
#' This construction is implemented in `add_linear_approx_Q_function()`.
#'
#' ## Helper Functions
#'
#' The following helper functions for using approximation are available:
#'
#' * `approx_Q_value()` calculates approximate Q values given the weights in
#'    the model or specified weights.
#' * `approx_greedy_action()` uses approximate Q values given the weights in
#'    the model or specified weights to find the the greedy action for a
#'    state.
#' * `approx_greedy_policy()` calculates the greedy-policy
#'    for the approximate Q values given the weights in
#'    the model or specified weights.
#'
#' ## Episodic Semi-gradient Sarsa
#'
#' The implementation follows the algorithm given in Sutton and Barto (2018).
#'
#' @family solver
#' @family MDPTF
#'
#' @references
#' Sutton, Richard S., and Andrew G. Barto. 2018. Reinforcement Learning: An Introduction. Second. The MIT Press. [http://incompleteideas.net/book/the-book-2nd.html](http://incompleteideas.net/book/the-book-2nd.html).
#'
#' Alborz Geramifard, Thomas J. Walsh, Stefanie Tellex, Girish Chowdhary, Nicholas Roy, and Jonathan P. How. 2013. A Tutorial on Linear Function Approximators for Dynamic Programming and Reinforcement Learning. Foundations and Trends in Machine Learning 6(4), December 2013, pp. 375-451. \doi{10.1561/2200000042}
#'
#' @examples
#' # Example 1: A maze without walls. The step cost is 1. The start is top-left and
#' # the goal (+100 reward) is bottom-right.
#' # This is the ideal problem for a linear approximation of the Q-function
#' # using the x/y location as state features.
#'
#' m <- gw_maze_MDP(c(5, 5), start = "s(1,1)", goal = "s(5,5)")
#'
#' # construct state features as the x/y coordinates in the gridworld
#' state_features <- gw_s2rc(S(m))
#' state_features
#'
#' m <- add_linear_approx_Q_function(m, state_features)
#'
#' # constructed state-action features (X) and approximate Q function
#' # and gradient
#' m$approx_Q_function
#'
#' sol <- solve_MDP_APPROX(m, horizon = 1000, n = 10,
#'                      alpha = 0.01, epsilon = .5)
#' gw_plot(sol)
#' sol <- solve_MDP_APPROX(sol, horizon = 1000, n = 100,
#'           alpha = 0.01, epsilon = 0.05, verbose = FALSE, continue = TRUE)
#' gw_plot(sol)
#'
#' ### TESTS
#' sol <- solve_MDP_APPROX(m, horizon = 1000, n = 1,
#'                      alpha = 0.01, epsilon = .8, verbose = 2)
#' Q_values(sol)
#' gw_plot(sol)
#'
#' gw_animate(m, method = "APPROX:semi", horizon = 1000, n = 10,
#'                      alpha = 0.01, epsilon = .8)
#'
#' ###
#'
#' gw_matrix(sol, what = "value")
#'
#' # learned weights and state values
#' sol$solution$w
#'
#' # extracting approximate Q-values
#' approx_greedy_action(sol, "s(4,5)")
#' approx_Q_value(sol, "s(4,5)", "down")
#' approx_Q_value(sol)
#'
#' # extracting a greedy policy using the approximate Q-values
#' approx_greedy_policy(sol)
#'
#' # Example 2: Stuart Russell's 3x4 Maze using approximation
#' # The wall and the -1 absorbing state make linear approximation
#' # using just the position more difficult.
#' data(Maze)
#' gw_plot(Maze)
#'
#' Maze_approx <- add_linear_approx_Q_function(Maze)
#' sol <- solve_MDP_APPROX(Maze_approx, horizon = 100, n = 100,
#'                      alpha = 0.01, epsilon = 0.3)
#' gw_plot(sol)
#'
#' # Example 3: Use order-2 Fourier basis for features
#' order <- 2
#' cs <- expand.grid(0:order, 0:order)
#'
#' # convert state features to Fourier basis features
#' x <- state2features(S(Maze))
#' x
#'
#' x <- sweep(x, MARGIN = 2,
#'            STATS = apply(x, MARGIN = 2, max), FUN = "/")
#' x <- apply(cs, MARGIN = 1,
#'            FUN = function(c) cos(pi * x %*% c))
#' rownames(x) <- S(Maze)
#' x
#'
#' Maze_approx <- add_linear_approx_Q_function(Maze, state_features = x)
#'
#' sol <- solve_MDP_APPROX(Maze_approx, horizon = 100, n = 100, alpha = 0.1, epsilon = .5)
#' sol <- solve_MDP_APPROX(sol, horizon = 100, n = 1000, alpha = 0.1, epsilon = .1, continue = TRUE)
#' gw_plot(sol)
#' 
#' # order-2 fourier features can be specified easier
#' Maze_approx <- add_linear_approx_Q_function(Maze, 
#'       transformation = transformation_fourier(min = c(0,0), max = c(3,4), order = 2))
#' sol <- solve_MDP_APPROX(Maze_approx, horizon = 1000, n = 200, alpha = 0.1, epsilon = .1)
#' gw_plot(sol)
#' 
#' @inheritParams solve_MDP
#' @param method string; one of the following solution methods: `'semi_gradient_sarsa'`
#' @param alpha step size.
#' @param epsilon used for the \eqn{\epsilon}-greedy behavior policies.
#' @param n number of episodes used for learning.
#' @param w an initial weight vector. By default a vector with 0s is used.
#'
#' @inherit solve_MDP return
#'
#' @export
solve_MDP_APPROX <-
  function(model,
           method = "semi_gradient_sarsa",
           horizon = NULL,
           discount = NULL,
           alpha = 0.01,
           epsilon = 0.2,
           n = 1000,
           w = NULL,
           ...,
           matrix = TRUE,
           continue = FALSE,
           progress = TRUE,
           verbose = FALSE) {
    .nodots(...)
    
    ### FIXME: It should work with MDPF with states
    if (!inherits(model, "MDPE"))
      stop("This model needs to be a MDP environment.")
    
    method <-
      match.arg(method, c("semi_gradient_sarsa"))
    
    q <- model$approx_Q_function
    if (is.null(q$w_init) || is.null(q$f) || is.null(q$gradient))
      stop(
        "Approx q-function f, the gradient, or w_init are missing in the q_function element in the model!"
      )
    
    # this code solve MDP and MDPTF
    model <- .prep_model(model,
                         horizon,
                         discount,
                         matrix,
                         verbose,
                         progress)
    
    if (verbose)
      progress <- FALSE
    
    
    if (progress) {
      pb <- my_progress_bar(n + 1L, name = "solve_MDP")
      pb$tick(0)
    }
    
    # Initialize w
    if (continue) {
      if (!is.null(w))
        stop("continue and w cannot be both specified!")
      w <- model$solution$w
    } else {
      w <- w %||% q$w_init
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
      
      model$solution <- list(
        method = method,
        alpha = alpha,
        epsilon = epsilon,
        n = n,
        w = w,
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
      cat("\nalpha:            ", deparse(alpha))
      cat("\nepsilon:          ", epsilon)
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
      
      # MDP: state is an id, MDPTF: state is a features
      s <- start(model)
      
      a <- approx_greedy_action(model, s, w, epsilon)
      
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
          
          w <- w + alpha * (r - q$f(s_features, a, w)) * q$gradient(s_features, a, w)
          
          if (verbose > 1) {
            print(w)
            cat(sprintf("New q(s,a)=%.3f\n", q$f(s_features, a, w)))
          }
          
          break
        }
        
        s_prime_features <- normalize_state_features(s_prime, model)
        a_prime <- approx_greedy_action(model, s_prime, w, epsilon)
        
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
          stop("Temporal differencing error becomes too large which indicates instability. Reduce alpha.")
        
        w <- w + alpha * td_error * q$gradient(s_features, a, w)
        
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

#' @rdname solve_MDP_APPROX
#' @export
add_linear_approx_Q_function <- function(model, ...) {
  UseMethod("add_linear_approx_Q_function")
}

#' @rdname solve_MDP_APPROX
#' @param state_features a matrix with state features. Each row is the feature
#'    vector for a state.
#' @param transformation a transformation function that is applied to the feature vector
#'    before it is used in the linear approximation.
#' @export
add_linear_approx_Q_function.MDP <- function(model,
                                             state_features = NULL,
                                             transformation = NULL,
                                             ...) {
  .nodots(...)
  state_features <- state_features %||% state2features(S(model))
  
  if (!is.matrix(state_features) ||
      !nrow(state_features) == c(length(S(model))))
    stop("state_features needs to be a matrix with one row per state!")
  
  model$state_features <- state_features
  
  # apply any transformation
  if (is.null(transformation))
    dim_s <- ncol(state_features)
  else
    dim_s <- length(transformation(state_features[1L, ]))
  
  n_A <- length(A(model))
  
  # for MDPTF we get state feature. Convert to action-state feature
  x <- function(s, a) {
    # use a transformation function
    if (!is.null(transformation))
      s <- transformation(s)
    
    a <- normalize_action_id(a, model)
    
    # 1 intercept
    #x <- numeric(1L + n_A * dim_s)
    #x[1L] <- 1
    #a_pos <- 2L + (a - 1L) * dim_s
    #x[a_pos:(a_pos + dim_s - 1L)] <- s
    
    # no intercept
    #x <- numeric(n_A * dim_s)
    #a_pos <- 1L + (a - 1L) * dim_s
    #x[a_pos:(a_pos + dim_s - 1L)] <- s
    
    # intercept per action
    x <- numeric(n_A * (1L + dim_s))
    a_pos <- 2L + (a - 1L) * (dim_s + 1L)
    x[a_pos - 1L] <- 1
    x[a_pos:(a_pos + dim_s - 1L)] <- s
    
    x
  }
  
  model$approx_Q_function <- list(
    f = function(s, a, w)
      sum(w * x(s, a)),
    gradient = function(s, a, w)
      x(s, a),
    #w_init = numeric(1L + n_A * dim_s)
    #w_init = numeric(n_A * dim_s)
    w_init = numeric(n_A * (1L + dim_s))
  )
  
  model
}


#' @rdname solve_MDP_APPROX
#' @export
add_linear_approx_Q_function.MDPTF <- function(model, transformation = NULL, ...) {
  .nodots(...)
  if (is.null(transformation))
    dim_s <- length(model$start)
  else
    dim_s <- length(transformation(model$start))
  
  n_A <- length(A(model))
  
  # for MDPTF we get state feature. Convert to action-state feature
  x <- function(s, a) {
    # use a transformation function
    if (!is.null(transformation))
      s <- transformation(s)
    
    a <- normalize_action_id(a, model)
    
    # 1 intercept
    #x <- numeric(1L + n_A * dim_s)
    #x[1L] <- 1
    #a_pos <- 2L + (a - 1L) * dim_s
    #x[a_pos:(a_pos + dim_s - 1L)] <- s
    
    # no intercept
    #x <- numeric(n_A * dim_s)
    #a_pos <- 1L + (a - 1L) * dim_s
    #x[a_pos:(a_pos + dim_s - 1L)] <- s
    
    # intercept per action
    x <- numeric(n_A * (1L + dim_s))
    a_pos <- 2L + (a - 1L) * (dim_s + 1L)
    x[a_pos - 1L] <- 1
    x[a_pos:(a_pos + dim_s - 1L)] <- s
    
    x
  }
  
  model$approx_Q_function <- list(
    f = function(s, a, w)
      sum(w * x(s, a)),
    gradient = function(s, a, w)
      x(s, a),
    #w_init = numeric(1L + n_A * dim_s)
    #w_init = numeric(n_A * dim_s)
    w_init = numeric(n_A * (1L + dim_s))
  )
  
  model
}

#' @rdname solve_MDP_APPROX
#' @param min,max vectors with the minimum and maximum values for each feature. 
#'    This is used to scale the feature to the \eqn{[0,1]} interval for the 
#'    Fourier basis.
#' @param order order for the Fourier basis.
#' @param cs an optional matrix or a data frame to specify cs values for the 
#'    Fourier features selectively (overrides `order`).
#' @export
transformation_fourier <- function(min, max, order, cs = expand.grid(0:order, 0:order)) {
  function(x) {
    x <- (x - min) / max
    apply(cs, MARGIN = 1,
          FUN = function(c) cos(pi * x %*% c))
  }
}
  
#' @rdname solve_MDP_APPROX
#' @param state a state (index or name)
#' @param action an action (index or name)
#' @param w a weight vector
#' @export
approx_Q_value <- function(model,
                           state = NULL,
                           action = NULL,
                           w = NULL) {
  action <- action %||% A(model)
  
  state <- state %||% S(model)
  if (is.null(state))
    stop("a state needs to be specified!")
  
  if (length(action) > 1L)
    return(sapply(
      action,
      FUN = function(a)
        approx_Q_value(model, state, a, w)
    ))
  
  w <- w %||% model$solution$w
  if (is.null(w))
    stop("weight vector w is missing.")
  
  sf <- normalize_state_features(state, model)
  if (nrow(sf) == 1L)
    model$approx_Q_function$f(drop(sf), action, w)
  else
    apply(sf, MARGIN = 1, model$approx_Q_function$f, action, w)
}

#' @rdname solve_MDP_APPROX
#' @export
approx_greedy_action <- function(model,
                                 state,
                                 w = NULL,
                                 epsilon = 0) {
  if (epsilon == 0 || runif(1) > epsilon) {
    a <- which.max.random(approx_Q_value(model, state, w  = w))
  } else {
    a <- sample.int(length(A(model)), 1L)
  }
  
  normalize_action(a, model)
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
