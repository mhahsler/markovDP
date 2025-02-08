#' Available Actions in a State
#'
#' Determine the set of actions available in a state.
#'
#' Unavailable actions are modeled as actions that have an immediate
#' reward of `-Inf` in the reward function. For a maze, also actions that
#' do not change the state can be considered unavailable.
#' @name available_actions
#' @family MDP
#'
#' @param model a [MDP] object.
#' @param state a character vector specifying the states.
#' @param neg_inf_reward logical; consider an action that produced `-Inf` reward to all
#'  end states unavailable?
#' @param stay_in_place logical; consider an action that results in the same state
#'  with a probability of 1 as unavailable. Note that this means that
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
#' available_actions(DynaMaze, state = "s(1,1)", stay_in_place = TRUE)
#'
#' # An action that leaves the grid currently is allowed but does not do
#' # anything.
#' act(DynaMaze, "s(1,1)", "up")
#' @returns a vector with the available actions.
#' @export
available_actions <- function(model,
                                  state,
                                  neg_inf_reward = TRUE,
                                  stay_in_place = FALSE) {
  UseMethod("available_actions")
}

#' @export
available_actions.MDP <- function(model,
                              state,
                              neg_inf_reward = TRUE,
                              stay_in_place = FALSE) {
  
  # calculate using the transition probabilities
  # deal with a single state
  if (neg_inf_reward) {
    if (is.null(state)) {
      acts_reward <- sapply(A(model), function(a) {
        rowSums(reward_matrix(model, a) != -Inf & transition_matrix(model, a) > 0) > 0
      })
      
    } else if (length(state) == 1L) {
      acts_reward <- sapply(A(model), function(a) {
        sum(reward_matrix(model, a, state) != -Inf & transition_matrix(model, a, state) > 0) > 0
      })
      acts_reward <- rbind(acts_reward)
      rownames(acts_reward) <- normalize_state(state, model)
        
    } else {
      acts_reward <- sapply(A(model), function(a) {
        rowSums(reward_matrix(model, a, state) != -Inf & transition_matrix(model, a, state) > 0) > 0
      })
    }
  }
 
  if (stay_in_place) {
    acts_stay <- t(sapply(state, FUN = function(s) 
      transition_matrix(
        model,
        action = A(model),
        start.state = s,
        end.state = s,
        simplify = TRUE) != 1
    ))
    rownames(acts_stay) <- normalize_state(state, model)
  }
  
  if (neg_inf_reward && stay_in_place)
    acts <- acts_reward & acts_stay
  else if (neg_inf_reward && !stay_in_place)
    acts <- acts_reward
  else if (!neg_inf_reward && stay_in_place)
    acts <- acts_stay
  else 
    stop ("No available action rule selected!")
  
  acts
}


#' @export
available_actions.MDPTF <- function(model,
                              state,
                              neg_inf_reward = TRUE,
                              stay_in_place = FALSE) {
  # figure them out by trying. This may not be perfect!
  A <- A(model)
  res <- lapply(A, FUN = function (a) act(model, state, a))
  
  if (neg_inf_reward) rew <- sapply(res, "[[", 1L) != -Inf
  else rew <- rep(TRUE, length(A))
  
  if (stay_in_place) sp <- sapply(res, "[[", 2L) != state
  else sp <- rep(TRUE, length(A))
  
  normalize_action(A[rew & sp], model) 
}

