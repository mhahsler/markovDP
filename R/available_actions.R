#' Available Actions in a State
#'
#' Determine the set of actions available in a state.
#'
#' Unavailable actions are modeled as actions that have an immediate
#'  reward of `-Inf` in the reward function.
#' @name available_actions
#' @family MDP
#'
#' @param model a [MDP] object.
#' @param state a character vector of length one specifying the state.
#' @param neg_inf_reward logical; consider an action that produced `-Inf` reward to all
#'  end states unavailable?
#' @param stay_in_place logical; consider an action that results in the same state
#'  with a probability of 1 as unavailable. Note that this will mean that
#'  absorbing states have no available action!
#' @returns a character vector with the available actions.
#'
#' @author Michael Hahsler
#' @examples
#' data(DynaMaze)
#' gw_plot(DynaMaze)
#'
#' # The the following actions are always available:
#' DynaMaze$actions
#' 
#' # only right and down is unavailable for s(1,1) because they
#' #   make the agent stay in place.
#' available_actions(DynaMaze, state = "s(1,1)")
#'
#' # An action that leaves the grid currently is allowed but does not do
#' # anything.
#' act(DynaMaze, "s(1,1)", "up")
#' @returns a vector with the available actions.
#' @export
available_actions <- function(model,
                              state,
                              neg_inf_reward = TRUE,
                              stay_in_place = TRUE) {
  acts <- model$actions
  
  if (stay_in_place) {
    acts <- acts[transition_matrix(
      model,
      action = acts,
      start.state = state,
      end.state = state,
      simplify = TRUE
    ) != 1L]
    
    # absorbing state 
    if (length(acts) < 1L)
      acts <- model$actions
  }
  
  if (neg_inf_reward) {
    acts <- acts[!apply(reward_matrix(model, action = acts, start.state = state, simplify = TRUE) == -Inf,
                        MARGIN = 1,
                        all)]
  }
  
  .normalize_action(acts, model)
}