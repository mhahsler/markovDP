# See True online Sarsa(lambda): Sutton & Barto 2018, p. 307

#' Eligibility Trace Algorithms with Linear Function Approximation
#'
#' Implementation of Sarsa(\eqn{\lambda}) and GTD(\eqn{\lambda}) following Sutton & Barto (2018).
#' 
#' Approximation functions are added to models using 
#' `add_linear_approx_Q_function()`.
#' 
#' ## Implemented methods
#' 
#' * `'Sarsa'` - True online Sarsa(\eqn{\lambda}).
#' * `'GTD'` - GTD(\eqn{\lambda}).
#' 
#' ## Schedules
#' 
#' * epsilon schedule: `t` is increased by each processed episode.
#' * alpha schedule: `t` is increased by each processed episode.
#' 
#' @family solver
#' @family MDP
#' @family MDPTF
#'
#' 
#' @inheritParams solve_MDP_APPROX
#' @param method string; one of the following solution methods:
#'    * `'Sarsa'` - True online Sarsa(\eqn{\lambda}).
#'    * `'GTD'` - GTD(\eqn{\lambda}).
#' @param lambda the trace-decay parameter for the an accumulating trace.
#'
#' @inherit solve_MDP return
#'
#' @author Michael Hahsler
#' @references
#' Sutton, Richard S., and Andrew G. Barto. 2018. Reinforcement Learning: An Introduction. Second. The MIT Press. [http://incompleteideas.net/book/the-book-2nd.html](http://incompleteideas.net/book/the-book-2nd.html).
#'
#' @examples
#' 
#' # Example 1: A Maze without walls
#' Maze <- gw_maze_MDP(c(5, 5), start = "s(1,1)", goal = "s(5,5)") 
#' 
#' Maze_approx <- add_linear_approx_Q_function(Maze)
#' sol <- solve_MDP_APPROX_LAMBDA(Maze_approx, horizon = 100, n = 100,
#'                      lambda = .1, 
#'                      alpha = schedule_exp(0.2, .1), 
#'                      epsilon = schedule_exp(1, .1)
#'                      )
#' gw_plot(sol)
#' gw_matrix(sol, what = "value")
#' approx_V_plot(sol)
#' 
#' # Example 2: Stuart Russell's 3x4 Maze
#' data(Maze)
#' Maze_approx <- add_linear_approx_Q_function(Maze, 
#'                         transformation = transformation_fourier_basis, 
#'                         order = 2)
#' sol <- solve_MDP_APPROX_LAMBDA(Maze_approx, horizon = 100, n = 200,
#'                      lambda = 0.1, 
#'                      alpha = schedule_exp(0.2, .1), 
#'                      epsilon = schedule_exp(1, .1))
#' gw_plot(sol)
#' gw_matrix(sol, what = "value")
#' approx_V_plot(sol)
#' 
#' @export
solve_MDP_APPROX_LAMBDA <-
  function(model,
           method = "sarsa",
           horizon = NULL,
           discount = NULL,
           lambda = 0.1,
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
        
        # get features x and x'
      
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

