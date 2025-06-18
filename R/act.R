#' Perform an Action
#'
#' Performs an action in a state and returns the new state and reward.
#'
#' @family MDP
#' @family MDPTF
#'
#' @param model an MDP model.
#' @param state the current state.
#' @param action the chosen action. If the action is not specified (`NULL`) and
#'   the MDP model contains a policy, then the action is chosen according to the
#'   policy.
#' @param fast logical; if `TRUE` then extra state id to label conversions are avoided.
#' @param ... if action is unspecified, then the additional parameters are
#'   passed on to [action()] to determine the action using the model's policy.
#'
#' @returns a names list with
#'  the `reward` and the next `state_prime`.
#'
#' @author Michael Hahsler
#' @examples
#' data(Maze)
#'
#' act(Maze, "s(1,3)", "right")
#'
#' # solve the maze and then ask for actions using the policy
#' sol <- solve_MDP(Maze)
#' act(sol, "s(1,3)")
#'
#' # make the policy in sol epsilon-soft and ask 10 times for the action
#' replicate(10, act(sol, "s(1,3)", epsilon = .2))
#' @export
act <- function(model, state, action, fast = FALSE, ...) {
  UseMethod("act")
}

#' @rdname act
#' @export
act.MDP <- function(model,
                    state,
                    action = NULL,
                    fast = FALSE,
                    ...) {
  # follow the policy in the model?
  action <- action %||% action(model, state, ...)
  
  #action <- normalize_action(action, model)
  #state <- normalize_state(state, model)
  
  sp <- sample.int(length(S(model)), 1L, prob = transition_matrix(model, action, state, sparse = FALSE))
  
  if (!fast)
    sp <- normalize_state(sp, model)
  
  r <-  reward_matrix(model, action, state, sp)
  
  list(reward =  r, state_prime = sp)
}

# fast is ignored since we used state features
#' @rdname act
#' @export
act.MDPTF <- function(model, state, action, fast = FALSE, ...) {
  state <- normalize_state_features(state, model)
  spr <- model$transition_func(model, state, action)

  if (!fast && !is.null(model$states))
    spr$state_prime <- normalize_state(spr$state_prime, model)
  
  spr
  }
