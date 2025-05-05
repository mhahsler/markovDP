#' @rdname policy_evaluation
#' @param ... further parameters are passed on to [`sample_MDP()`].
#' @param n number of simulated episodes.
#' @param horizon maximum horizon for episodes.
#' @param first_visit if `TRUE` then only the first visit of a state/action pair
#'   in an episode is used to update Q, otherwise, every-visit update is used.
#' @export
policy_evaluation_MC <- function(model,
           pi = NULL,
           n = 1000,
           horizon = NULL,
           first_visit = TRUE,
           ...,
           progress = TRUE,
           verbose = FALSE) {

  S <- S(model)
  A <- A(model)
  horizon <- horizon %||% model$horizon
  discount <- model$discount
  
  if(!is.null(pi))
    model <- add_policy(model, pi)
  
  V <- V_zero(model)
  V_N <- integer(length(V))
  
  if (progress)
    pb <- my_progress_bar(n + 1L, name = "solve_MDP")
  
  on.exit({
    if (progress) {
      pb$tick(0)
      pb$terminate()
    }
    
    if (e < n)
      warning("Manual interupt: MDP solver stopped at episode ", e)
    
    if (verbose) {
      cat("\nTerminated after episode:", e, "\n")
    }
    
    V <- V/V_N
    V[is.na(V)] <- 0
    return(V)
  })
  
  # Loop through N episodes
  e <- 0L
  while (e < n) {
    e <- e + 1L
    if (progress)
      pb$tick()
    
    # follow the policy!
    ep <- sample_MDP(
      model,
      n = 1,
      horizon = horizon,
      epsilon = 0,
      exploring_starts = FALSE,
      trajectories = TRUE,
      progress = FALSE,
      verbose = FALSE,
      ...
    )$trajectories
    
    if (verbose > 1) {
      cat(paste(
        "\n****************** Episode",
        e,
        "******************\n"
      ))
      print(ep)
      cat("\n")
    }
    
    G <- 0
    for (i in rev(seq_len(nrow(ep)))) {
      r_t_plus_1 <- ep$r[i]
      s_t <- normalize_state_id(ep$s[i], model)
      a_t <- ep$a[i]
      
      G <- discount * G + r_t_plus_1
      
      # Only update for first visit of a s
      if (first_visit && i < 2L &&
          (any(s_t == ep$s[1:(i - 1L)])))
        next
      
      # V <- avg(Returns)
      # running average instead of averaging Returns lists.
      V_N[s_t] <-  V_N[s_t] + 1L
      V[s_t] <- V[s_t] + G
    }
    
  } 
    # return is handled by on.exit()
}
  