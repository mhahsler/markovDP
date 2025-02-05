

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
#' @family MDPE
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
#' m <- gw_random_maze(5, wall_prob = 0)
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
#'                      alpha = 0.01, epsilon = .7)
#' sol <- solve_MDP_APPROX(sol, horizon = 1000, n = 100, 
#'           alpha = 0.01, epsilon = 0.1, verbose = FALSE, continue = TRUE)
#' 
#' gw_plot(sol)
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
#'
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
#' sol <- solve_MDP_APPROX(sol, horizon = 100, n = 100, alpha = 0.1, epsilon = .1, continue = TRUE)
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
    
    method <-
      match.arg(method, c("semi_gradient_sarsa"))
    
    q <- model$approx_Q_function
    if (is.null(q$w_init) || is.null(q$f) || is.null(q$gradient))
      stop(
        "Approx q-function f, the gradient, or w_init are missing in the q_function element in the model!"
      )
    
    # this code solve MDP and MDPE
    model <- .prep_model(model,
                         horizon,
                         discount,
                         matrix,
                         verbose,
                         progress,
                         allow_MDPE = TRUE)
   
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
        policy = list(approx_greedy_policy(model, w))
      )
      return(model)
    })
    
    if (verbose) {
      cat("Running", method)
      cat("\nalpha:            ", deparse(alpha))
      cat("\nepsilon:          ", epsilon)
      cat("\nn                 ", n, "\n")
    }
    
    A <- A(model)
    discount <- model$discount
    
    # start: sample only if we have more than one start state
    start_sample <- FALSE
    if (inherits(model, "MDPE")) {
      start <- model$start
    } else { # MDP
      start <- start_vector(model, sparse = "index")
      if (length(start) > 1L) {
        S <- S(model)
        start_prob <- start_vector(model, sparse = FALSE)
        start_sample <- TRUE
      }
    }
   
    # Note: for MDP s is a state id, but for MDPE s is a feature vector!
     
    # loop episodes
    e <- 0L
    while (e < n) {
      e <- e + 1L
      if (progress)
        pb$tick()
      
      if (start_sample) {
        s <- sample.int(length(S), 1L, prob = start_prob)
      } else {
        s <- start
      }
      
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
        
        # reached terminal state?
        if (absorbing_states(model, state = s_prime)) {
          if (verbose > 1) {
            cat(
              sprintf(
                "\nStep %i - s=%s a=%i r=%.3f s'=%s a'=%i q(s,a)=%.3f error=%.3f\n",
                i,
                s,
                a,
                r,
                s_prime,
                a_prime,
                q$f(s, a, w),
                r - q$f(s, a, w)
              )
            )
            print(w)
          }
          
          w <- w + alpha * (r - q$f(s, a, w)) * q$gradient(s, a, w)
          
          if (verbose > 1) {
            print(w)
            cat(sprintf("New q(s,a)=%.3f\n", q$f(s, a, w)))
            cat("Terminal state reached!")
          }
          
          break
        }
        
        a_prime <- approx_greedy_action(model, s_prime, w, epsilon)
        
        if (verbose > 1) {
          if (i == 1L) {
            cat("\n*** Episode", e, "***\n")
          }
          cat(
            sprintf(
              "\nStep %i - s=%s a=%i r=%.3f s'=%s a'=%i  q(s,a)=%.3f error=%.3f\n",
              i,
              s,
              a,
              r,
              s_prime,
              a_prime,
              q$f(s, a, w),
              r + discount * q$f(s_prime, a_prime, w) - q$f(s, a, w)
            )
          )
          print(w)
        }
        
        # this is one-step Sarsa
        # unavailable actions return a reward of -r
        if (!is.finite(r))
          stop("Inf reward not supported!")
        
        w <- w + alpha * (r + discount * q$f(s_prime, a_prime, w) - q$f(s, a, w)) * q$gradient(s, a, w)
        
        if (verbose > 1) {
          print(w)
          cat(sprintf("New q(s,a)=%.3f\n", q$f(s, a, w)))
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
add_linear_approx_Q_function.MDP <- function(model, state_features = NULL, transformation = NULL, ...) {
  .nodots(...)
  state_features <- state_features %||% state2features(S(model))
  
  if (!is.matrix(state_features) ||
      !nrow(state_features) == c(length(S(model))))
    stop("state_features needs to be a matrix with one row per state!")
  
  # apply any transformation
  if (!is.null(transformation))
    state_features <- t(apply(state_features, MARGIN = 1, transformation))
  
  A <- A(model)
  
  # Construct state-action feature vectors as a list of actions
  X <- lapply(
    seq_along(A),
    FUN = function(a) {
      x <- matrix(
        0,
        nrow = nrow(state_features),
        ncol = 1L + length(A) * ncol(state_features),
        dimnames = list(rownames(state_features), NULL)
      )
      x[, 1L] <- 1
      a_pos <- 2L + (a - 1L) * ncol(state_features)
      x[, a_pos:(a_pos + ncol(state_features) - 1L)] <- state_features
      x
    }
  )
  names(X) <- A
  
  # The state-action value is approximated by the inner product t(w)x
  f <- function(s, a, w)
    sum(w * X[[a]][s, ])
  
  # The gradient for linear approximation is just x
  gradient <- function(s, a, w)
    X[[a]][s, ]
  
  # Initialize a new weight vector. We use a set of weights for x and y
  # for each action.
  w_init <- numeric(ncol(X[[1]]))
  
  model$approx_Q_function <- list(
    f = f,
    gradient = gradient,
    w_init = w_init,
    X = X
  )
  model
}


#' @rdname solve_MDP_APPROX
#' @export
add_linear_approx_Q_function.MDPE <- function(model, transformation = NULL, ...) {
  .nodots(...)
  if (is.null(transformation))
    dim_s <- length(model$start)
  else
    dim_s <- length(transformation(model$start))
  
  n_A <- length(A(model))
  
  # state feature to action-state feature
  x <- function(s, a) {
    #cat("s:")
    #print(s)
    #cat("a:")
    #print(a)
    
    # use a transformation function
    if (!is.null(transformation))
      s <- transformation(s)
    
    a <- normalize_action_id(a, model)
    x <- numeric(1L + n_A * dim_s)
    x[1L] <- 1
    a_pos <- 2L + (a - 1L) * dim_s
    x[a_pos:(a_pos + dim_s - 1L)] <- s
    
    #cat("x:")
    #print(x)
    #cat("\n\n")
    
    x
  }
  
  model$approx_Q_function <- list(
    f = function(s, a, w) sum(w * x(s, a)), 
    gradient = function(s, a, w) x(s, a), 
    w_init = numeric(1L + n_A * dim_s) )
  model
}

#' @rdname solve_MDP_APPROX
#' @param state a state (index or name)
#' @param action an action (index or name)
#' @param w a weight vector
#' @export
approx_Q_value <- function(model,
                           state,
                           action = NULL,
                           w = NULL) {
 
 ### whoe table comes from q_values!!!!
   
  action <- action %||% A(model)
  if (length(action) > 1L)
    return(sapply(
      action,
      FUN = function(a)
        approx_Q_value(model, state, a, w)
    ))
  
  w <- w %||% model$solution$w
  if (is.null(w))
    stop("weight vector w is missing.")
  
  state <- normalize_state_label(state, model)
  model$approx_Q_function$f(state2features(state), action, w)
}

#' @rdname solve_MDP_APPROX
#' @export
approx_greedy_action <- function(model,
                                 state,
                                 w = NULL,
                                 epsilon = 0) {
  state <- normalize_state_label(state, model)
  
  # calculate q(s, a) for all a and pick the largest
  w <- w %||% model$solution$w
  if (is.null(w))
    stop("weight vector w is missing.")
  
  q <- model$approx_Q_function
  qs <- sapply(
    seq_along(A(model)),
    FUN = function(a)
      q$f(state2features(state), a, w)
  )
  
  if (epsilon == 0 || runif(1) > epsilon) {
    a <- which.max.random(qs)
  } else {
    a <- sample.int(length(A(model)), 1L)
  }
  
  normalize_action(a, model)
}

#' @rdname solve_MDP_APPROX
#' @export
approx_greedy_policy <- function(model, w = NULL) {
  w <- w %||% model$solution$w
  acts <- sapply(
    S(model),
    FUN = function(s)
      approx_greedy_action(model, s, w, epsilon = 0)
  )
  
  data.frame(
    state = S(model),
    V = sapply(
      seq_along(S(model)),
      FUN = function(i)
        approx_Q_value(model, i, acts[i], w)
    ),
    action = acts,
    row.names = NULL
  )
}
