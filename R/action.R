#' Choose an Action Given a Policy
#'
#' Returns an action given a deterministic policy.
#' The policy can be made epsilon-soft.
#'
#' @family policy

#' @param model a solved [MDP].
#' @param state the state.
#' @param epoch what epoch of the policy should be used. Use 1 for converged policies.
#' @param epsilon make the policy epsilon soft.
#' @param ... further parameters are passed on.
#' @return The name of the optimal action as a factor.
#' @author Michael Hahsler
#' @examples
#' data("Maze")
#' Maze
#'
#' sol <- solve_MDP(Maze)
#' policy(sol)
#'
#' action(sol, state = "s(1,3)")
#'
#' ## choose from an epsilon-soft policy
#' table(replicate(100, action(sol, state = "s(1,3)", epsilon = 0.1)))
#' @export
action <-
  function(model,
           state,
           epsilon = 0,
           epoch = 1,
           ...) {
    ### TODO: FIXME - Available actions only?
    A <- A(model)
    
    if (epsilon == 1)
      return(normalize_action(sample.int(length(A), size = 1L), model))
    
    
    if (is.matrix(state) && nrow(state) != 1L ||
        is.vector(state) && length(state) != 1)
      stop("action() requires a single state!")
    
    # we can have a policy for tabular solutions or
    #     an approximation (where we choose a greedy action).
    if (!is.null(model$solution$policy[[1]])) {
      a <- policy(model, epoch)$action[normalize_state_id(state, model)]
    } else if (!is.null(model$solution$w)) {
      a <- approx_greedy_action(model, state)
    } else
      stop("No explicit policy or approximate Q-function available! Set epsilon to 1 for random actions.")
    
    if (epsilon == 0)
      return(a)
    
    if (length(A) > 1L && runif(1) < epsilon) {
      a <- sample.int(length(A), size = 1L)
      a <- normalize_action(a, model)
    }
    
    a
  }
