
#' @rdname policy_evaluation
#' @param V a vector with estimated state values representing a value function.
#'    If `model` is a solved model, then the state
#'    values are taken from the solution.
#' @param pi a policy as a data.frame with at least columns for states and action. If `NULL`,
#'     then the policy in model is used.
#' @param k_backups number of look ahead steps used for approximate policy evaluation
#'    used by the policy iteration method. Set k_backups to `Inf` to only use
#'    \eqn{\theta} as the stopping criterion.
#' @param theta stop when the largest state Bellman error (\eqn{\delta = V_{k+1} - V}) 
#'    is less than \eqn{\theta}.
#' @export
policy_evaluation_bellman <-
  function(model,
           pi = NULL,
           V = NULL,
           k_backups = 1000L,
           theta = 1e-3,
           progress = TRUE,
           verbose = FALSE) {
    if (!inherits(model, "MDP")) {
      stop("'model' needs to be of class 'MDP'.")
    }
    
    if (is.null(pi))
      pi <- policy(model)
    
    if (progress) {
      pb <- my_progress_spinner(name = "policy_evaluation")
      pb$tick(0)
    }
    
    S <- S(model)
    A <- A(model)
    gamma <- model$discount %||% 1
    
    # start with all 0s if no previous U is given
    if (is.null(V)) {
      V <- V_zero(model)
    }
    
    # we cannot count more than integer.max
    if (k_backups > .Machine$integer.max) {
      k_backups <- .Machine$integer.max
      warning("Using the maximum number of backups of", k_backups)
    }
    
    # use r_pi and p_pi to make it faster and not do a complete Bellman update
    p_pi <- induced_transition_matrix(model, pi)
    r_pi <- induced_reward_matrix(model, pi)
    r_pi <- rowSums(p_pi * r_pi)
    
    for (i in seq_len(k_backups)) {
      if (progress)
        pb$tick()
      
      v <- V
      
      #V <- bellman_operator(model, pi, V)
      # do this faster directly
      V <- r_pi + gamma * p_pi %*% V
      
      delta <- max(abs(v - V), na.rm = TRUE)
      
      if (verbose) {
        cat("Backup step", i, ": delta =", delta, "\n")
      }
      
      if (delta < theta) {
        break
      }
    }
    
    if (progress)
      pb$terminate()
    
    V <- drop(V)
    names(V) <- S(model)
    V
  }




