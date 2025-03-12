# Explanation for how to parametrize q(s,a):
# https://danieltakeshi.github.io/2016/10/31/going-deeper-into-reinforcement-learning-understanding-q-learning-and-linear-function-approximation/

#' Solve MDPs with Temporal Differencing and Linear Function Approximation
#'
#' Solve the MDP control problem using state-value approximation
#' by semi-gradient Sarsa (temporal differencing) for
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
#' Linear approximation has a single optimum and can be optimized using
#' a simple update rule following the gradient of the state-action function
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
#' state-action feature vector is \eqn{\phi(s,a) = (0,0,0,1,3,4,0,0,0)}.
#' Each action component has three entries and the 1 represent the intercept
#' for the state feature vector. The zeros represent the components for the
#' two not chosen actions.
#'
#' The construction of the state-action values is implemented in `add_linear_approx_Q_function()`.
#'
#' The state feature function \eqn{\phi()} starts with raw state feature vector
#' \eqn{\mathbf{x} = (x_1,x_2, ..., x_m)} that
#' are either user-specified or constructed by parsing the state labels of
#' form `s(feature list)`.  Then an optional nonlinear transformation
#' can be performed (see [transformation]). 
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
#'
#' Alborz Geramifard, Thomas J. Walsh, Stefanie Tellex, Girish Chowdhary, Nicholas Roy, and Jonathan P. How. 2013. A Tutorial on Linear Function Approximators for Dynamic Programming and Reinforcement Learning. Foundations and Trends in Machine Learning 6(4), December 2013, pp. 375-451. \doi{10.1561/2200000042}
#'
#' Konidaris, G., Osentoski, S., & Thomas, P. 2011. Value Function Approximation in Reinforcement Learning Using the Fourier Basis. Proceedings of the AAAI Conference on Artificial Intelligence, 25(1), 380-385. \doi{10.1609/aaai.v25i1.7903}
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
#' set.seed(1000)
#' sol <- solve_MDP_APPROX(m, horizon = 1000, n = 100)
#' 
#' gw_plot(sol)
#' gw_matrix(sol, what = "value")
#'
#' # learned weights and state values
#' sol$solution$w
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
#' Maze_approx <- add_linear_approx_Q_function(Maze)
#' 
#' set.seed(1000)
#' sol <- solve_MDP_APPROX(Maze_approx, horizon = 100, n = 100,
#'                         alpha = schedule_exp(0.3, 0.01),
#'                         epsilon = schedule_exp(1, 0.1))
#' gw_plot(sol)
#' gw_matrix(sol, what = "value")
#' approx_V_plot(sol, res = 20)
#'
#'
#' # Example 3: Stuart Russell's 3x4 Maze using
#' #            order-1 Fourier basis for approximation
#'
#' Maze_approx <- add_linear_approx_Q_function(Maze,
#'       transformation = transformation_fourier_basis, order = 1)
#'
#' set.seed(1000)
#' sol <- solve_MDP_APPROX(Maze_approx, horizon = 100, n = 100,
#'                     alpha = schedule_exp(0.3, .01),
#'                     epsilon = schedule_exp(1, .1))
#'                     
#' gw_plot(sol)
#' gw_matrix(sol, what = "value")
#' approx_V_plot(sol, res = 20)
#'
#' # Example 4: Stuart Russell's 3x4 Maze using
#' #            order-2 Fourier basis for approximation
#' #            and eligibility traces: Sarsa(lambda)
#'
#' Maze_approx <- add_linear_approx_Q_function(Maze,
#'       transformation = transformation_fourier_basis, order = 1)
#'
#' set.seed(1000)
#' sol <- solve_MDP_APPROX(Maze_approx, horizon = 100, n = 100,
#'                     alpha = schedule_exp(0.3, .01),
#'                     epsilon = schedule_exp(1, .1),
#'                     lambda = 0.1)
#'                     
#' gw_plot(sol)
#' gw_matrix(sol, what = "value")
#' approx_V_plot(sol, res = 20)
#' 
#'  
#' @inheritParams solve_MDP_TD
#' @param method string; one of the following solution methods: `'sarsa'`
#' @param lambda the trace-decay parameter for the an accumulating trace. If `lambda = 0`
#'      then 1-step Sarsa is used.
#' @param w an initial weight vector. By default a vector with 0s is used.
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
           w = NULL,
           ...,
           matrix = TRUE,
           continue = FALSE,
           progress = TRUE,
           verbose = FALSE) {
    .nodots(...)
    if (lambda == 0)
      solve_MDP_APPROX_1_step(model,
                              method,
                              horizon,
                              discount,
                              alpha,
                              epsilon,
                              n,
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
           w = NULL,
           ...,
           matrix = TRUE,
           continue = FALSE,
           progress = TRUE,
           verbose = FALSE) {
    .nodots(...)
    
    if (!inherits(model, "MDPE"))
      stop("This model needs to be a MDP environment.")
    
    method <-
      match.arg(method, c("sarsa"))
    
    ### alpha/epsilon func
    alpha_seq <- NULL
    if (is.function(alpha))
      alpha_seq <- alpha(seq(1, n))
    else if(length(alpha) == n)
      alpha_seq <- alpha
    
    epsilon_seq <- NULL
    if (is.function(epsilon))
      epsilon_seq <- epsilon(seq(1, n))
    else if(length(epsilon) == n)
      epsilon_seq <- epsilon
    ###
    
    q <- model$approx_Q_function
    if (is.null(q$w_init) || is.null(q$f) || is.null(q$gradient))
      stop(
        "Approx q-function f, the gradient, or w_init are missing in the q_function element in the model!"
      )
    
    # this code solve MDP and MDPTF
    model <- .prep_model(model, horizon, discount, matrix, verbose, progress)
    
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
      
      ### alpha/epsilon func
      if (!is.null(epsilon_seq)) 
        epsilon <- epsilon_seq[e]
      if (!is.null(alpha_seq)) 
        alpha <- alpha_seq[e]
      ###
      
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
          stop(
            "Temporal differencing error becomes too large which indicates instability. Reduce alpha."
          )
        
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


solve_MDP_APPROX_lambda <-
  function(model,
           method = "sarsa",
           horizon = NULL,
           discount = NULL,
           alpha = schedule_exp(0.2, .1),
           epsilon = schedule_exp(1, .1),
           lambda = 0.1,
           n,
           w = NULL,
           ...,
           matrix = TRUE,
           continue = FALSE,
           progress = TRUE,
           verbose = FALSE) {
    .nodots(...)
    
    if (!inherits(model, "MDPE"))
      stop("This model needs to be a MDP environment.")
    
    methods <- c("Sarsa", "GTD")
    method <- match.arg(tolower(method), tolower(methods))
    method_id <- match(tolower(method), tolower(methods))
    
    ### alpha/epsilon func
    alpha_seq <- NULL
    if (is.function(alpha))
      alpha_seq <- alpha(seq(1, n))
    else if(length(alpha) == n)
      alpha_seq <- alpha
    
    epsilon_seq <- NULL
    if (is.function(epsilon))
      epsilon_seq <- epsilon(seq(1, n))
    else if(length(epsilon) == n)
      epsilon_seq <- epsilon
    ###
    
    
    qf <- model$approx_Q_function
    if (is.null(qf$w_init) || is.null(qf$f) || is.null(qf$gradient))
      stop(
        "Approx q-function f, the gradient, or w_init are missing in the q_function element in the model!"
      )
    
    # this code solve MDP and MDPTF
    model <- .prep_model(model, horizon, discount, matrix, verbose, progress)
    
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
      w <- w %||% qf$w_init
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
        lambda = lambda,
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
      
      
      ### alpha/epsilon func
      if (!is.null(epsilon_seq)) 
        epsilon <- epsilon_seq[e]
      if (!is.null(alpha_seq)) 
        alpha <- alpha_seq[e]
      ###
      
      # Initialize s and choose first action
      # MDP: state is an id, MDPTF: state is a features
      s <- start(model)
      a <- approx_greedy_action(model, s, w, epsilon)
      
      x <- qf$x(normalize_state_features(s, model), a)
      z <- qf$w_init 
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
        a_prime <- approx_greedy_action(model, s_prime, w, epsilon)
        
        # get x'
        x_prime <- qf$x(normalize_state_features(s_prime, model), a_prime)
        
        q <- sum(x * w)
        q_prime <- sum(x_prime * w)
        delta <- r + gamma * q_prime - q
        
        if (delta > 1e20)
          stop(
            "Temporal differencing error becomes too large which indicates instability. Reduce alpha or lambda."
          )
        
        z <- gamma * lambda * z + (1 - alpha + gamma * lambda * sum(z*x)) * x
        
        # Sarsa(lambda) .. on-policy true online
        if (method_id == 1L) {
          w <- w + alpha * (delta + q - q_old) * z - alpha * (q - q_old) * x
        }
        
        # GTD(lambda)
        else if (method_id == 2) {
          w <- w + alpha * delta * z - alpha * gamma * (1 - lambda) * sum(z*q) * x
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
          # print(w)
          print(z)
        }
        
        if (i >= horizon)
          break
        
        if(absorbing_states(Maze, state = s_prime))
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
#' @export
add_linear_approx_Q_function <- function(model, ...) {
  UseMethod("add_linear_approx_Q_function")
}

#' @rdname solve_MDP_APPROX
#' @param state_features a matrix with state features. Each row is the feature
#'    vector for a state.
#' @param transformation a transformation function that is applied to the feature vector
#'    before it is used in the linear approximation. See [transformation].
#' @export
add_linear_approx_Q_function.MDP <- function(model,
                                             state_features = NULL,
                                             transformation = transformation_linear_basis,
                                             ...) {
  state_features <- state_features %||% get_state_features(model)
  
  if (!is.matrix(state_features) ||
      !nrow(state_features) == c(length(S(model))))
    stop("state_features needs to be a matrix with one row per state!")
  
  if (is.null(colnames(state_features)))
    colnames(state_features) <- paste0("x", seq_len(ncol(state_features)))
  
  model$state_features <- state_features
  
  transformation <- transformation(model, ...)
  
  
  state_template <- transformation(state_features[1L, ])
  dim_s <- length(state_template)
  if (is.null(names(state_template)))
    names(state_template) <- paste0("x", seq_len(dim_s))
  n_A <- length(A(model))
  
  # s are state features. Convert to action-state feature.
  x <- function(s, a) {
    # use a transformation function
    s <- transformation(s)
    a <- normalize_action_id(a, model)
    
    # one component per action
    x <- numeric(n_A * dim_s)
    
    a_pos <- 1L + (a - 1L) * dim_s
    x[a_pos:(a_pos + dim_s - 1L)] <- s
    
    x
  }
  
  model$approx_Q_function <- list(
    x = x,
    f = function(s, a, w)
      sum(w * x(s, a)),
    gradient = function(s, a, w)
      x(s, a),
    #transformation = c(as.list(environment(transformation)), list(f = transformation)),
    transformation = transformation,
    w_init = setNames(numeric(n_A * dim_s), paste(
      rep(A(model), each = dim_s), names(state_template), sep = "."
    ))
  )
  
  model
}


#' @rdname solve_MDP_APPROX
#' @export
add_linear_approx_Q_function.MDPTF <- function(model,
                                               state_features = NULL,
                                               transformation = transformation_linear_basis,
                                               ...) {
  # MDPTF may have states
  if (!is.null(state_features))
    stop(
      "state_features cannot be specified for a MDPTF! The states are already in a factored representation."
    )
  
  # add state_features if we have a defined state space
  if (!is.null(S(model)))
    model$state_features <- state2features(S(model))
  
  # MDPTF has a start state matrix.
  # This is used to get min, max, etc for transformation
  start <- model$start[1, , drop = TRUE]
  if (is.null(names(start)))
    names(start) < paste0("x", seq_along(start))
  
  transformation <- transformation(model, ...)
  
  state_template <- transformation(start)
  dim_s <- length(state_template)
  if (is.null(names(state_template)))
    names(state_template) <- paste0("x", seq_len(dim_s))
  n_A <- length(A(model))
  
  # s are state features. Convert to action-state feature.
  x <- function(s, a) {
    # use a transformation function
    s <- transformation(s)
    a <- normalize_action_id(a, model)
    
    # one component per action
    x <- numeric(n_A * dim_s)
    
    a_pos <- 1L + (a - 1L) * dim_s
    x[a_pos:(a_pos + dim_s - 1L)] <- s
    
    x
  }
  
  model$approx_Q_function <- list(
    #x = x,
    f = function(s, a, w)
      sum(w * x(s, a)),
    gradient = function(s, a, w)
      x(s, a),
    transformation = transformation,
    w_init = setNames(numeric(n_A * dim_s), paste(
      rep(A(model), each = dim_s), names(state_template), sep = "."
    ))
  )
  
  model
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
  
  w <- w %||% model$solution$w
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
    model$approx_Q_function$f(drop(state), action, w)
  else
    apply(state, MARGIN = 1, model$approx_Q_function$f, action, w)
}

#' @rdname solve_MDP_APPROX
#' @export
approx_greedy_action <- function(model,
                                 state,
                                 w = NULL,
                                 epsilon = 0) {
  if (epsilon == 0 || runif(1) > epsilon) {
    a <- which.max.random(approx_Q_value(model, state, w = w))
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
