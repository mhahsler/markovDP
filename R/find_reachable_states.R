#' Find Reachable State Space from a Transition Model Function
#'
#' Finds the reachable state space from a transition model function
#' that takes the arguments `model`, `action`, `start.state` and
#' returns a named vector with the probabilities for the resulting end.states 
#' (typically only the ones with a probability greater than 1).
#'
#' The function performs a (depth-limited) depth-first traversal of the
#' search state space and returns a vector with the names of all encountered
#' states. This vector can be used as the states for creating a MDP model.
#'
#' @param transition_function a transition function (see details for requirements).
#' @param start_state labels of the start states.
#' @param actions a vector with the available actions.
#' @param model if needed, the model passed on to the transition model function. 
#' @param horizon only return states reachable in the given horizon.
#' @param progress logical; show a progress bar?
#' 
#' @returns a character vector with all reachable states.
#'
#' @author Michael Hahsler
#' 
#' @examples
#' # define a MDP for Tic-Tac-Toe
#' 
#' # state description: matrix with the characters _, x, and o
#' #                    can be converted into a label of 9 characters
#' 
#' # set of actions
#' A <- as.character(1:9)
#' 
#' # helper functions
#' ttt_empty_board <- function() matrix('_', ncol = 3, nrow = 3)
#' 
#' ttt_state2label <- function(state) paste(state, collapse = '')
#' 
#' ttt_label2state <- function(label) matrix(strsplit(label, "")[[1]], 
#'                                           nrow = 3, ncol = 3)
#' 
#' ttt_available_actions <- function(state) {
#'   if (length(state) == 1L) state <- ttt_label2state(state)
#'   which(state == "_")
#' }
#' 
#' ttt_result <- function(state, player, action) {
#'   if (length(state) == 1L) state <- ttt_label2state(state)
#'   
#'   if (state[action] != "_")
#'     stop("Illegal action.")
#'   state[action] <- player
#'   state
#' }
#' 
#' ttt_terminal <- function(state) {
#'   if (length(state) == 1L) state <- ttt_label2state(state)
#'   # Check the board for a win and return one of 
#'   # 'x', 'o', 'd' (draw), or 'n' (for next move)
#'   win_possibilities <- rbind(state, 
#'                              t(state), 
#'                              diag(state), 
#'                              diag(t(state)))
#' 
#'   wins <- apply(win_possibilities, MARGIN = 1, FUN = function(x) {
#'     if (x[1] != '_' && length(unique(x)) == 1) x[1]
#'     else '_'
#'   })
#' 
#'   if (any(wins == 'x')) 
#'     return('x')
#' 
#'   if (any(wins == 'o')) 
#'     return('o')
#' 
#'   # Check for draw
#'   if (sum(state == '_') < 1)
#'     return('d')
#' 
#'   return('n')
#' }
#' 
#' ttt_other <- function(player) if (player == 'x') 'o' else 'x'
#' 
#' # define the transition function
#' T <- function(model, action, start.state) {
#'   action <- as.integer(action)
#'   
#'   # absorbing states
#'   if (start.state %in% c('win', 'loss', 'draw')) {
#'     return(structure(1, names = start.state))
#'   }
#'   
#'   # avoid illegal action by making them a loss
#'   if (!(action %in% ttt_available_actions(ttt_label2state(start.state)))) {
#'     return(structure(1, names = "loss"))
#'   }
#'   
#'   # make x's move
#'   next_state <- ttt_result(ttt_label2state(start.state), 'x', action)
#'   
#'   # ttt_terminal?
#'   term <- ttt_terminal(next_state)
#'   
#'   if (term == 'x') {
#'     return(structure(1, names = "win"))
#'   } else if (term == 'o') {
#'     return(structure(1, names = "loss"))
#'   } else if (term == 'd') {
#'     return(structure(1, names = "draw"))
#'   }
#'   
#'   # it is o's turn
#'   actions_of_o <- ttt_available_actions(next_state)
#'   possible_end_states <- lapply(
#'     actions_of_o,
#'     FUN = function(a)
#'       ttt_result(next_state, 'o', a)
#'   )
#'   
#'   # fix ttt_terminal states
#'   term <- sapply(possible_end_states, ttt_terminal)
#'   possible_end_states <- sapply(possible_end_states, ttt_state2label)
#'   possible_end_states[term == 'x'] <- 'win'
#'   possible_end_states[term == 'o'] <- 'loss'
#'   possible_end_states[term == 'd'] <- 'draw'
#'   
#'   possible_end_states <- unique(possible_end_states)
#'   
#'   return(structure(rep(
#'     1 / length(possible_end_states), length(possible_end_states)
#'   ), names = possible_end_states))
#' }
#' 
#' # define the reward
#' R <- rbind(
#'   R_(                    value = 0),
#'   R_(end.state = 'win',  value = +1),
#'   R_(end.state = 'loss', value = -1),
#'   R_(end.state = 'draw', value = +.5),
#'   # Note: there is no more reward once the agent is in a terminal state
#'   R_(start.state = 'win',  value = 0),
#'   R_(start.state = 'loss', value = 0),
#'   R_(start.state = 'draw', value = 0)
#' )
#' 
#' # start state
#' start <- "_________"
#' 
#' # find the reachable state space
#' S <- find_reachable_states(T, start_state = start, actions = A) 
#' head(S)
#' 
#' tictactoe <- MDP(S, A, T, R, discount = 1, start = start, name = "TicTacToe")
#' tictactoe
#' @export
find_reachable_states <- function(transition_function, 
                                  start_state, 
                                  actions,
                                  model = NULL,
                                  horizon = Inf,
                                  progress = TRUE
                                  ) {
  states <- new.env(hash = TRUE)
  depth <- 0L
  
  if (is.null(model)) {
    model <- list(actions = actions, 
                  start = start_state, 
                  transition_prob = transition_function)
  }
  
  if (progress) {
    pb <- my_progress_spinner(name = "find_reachable_states", ticks = "reached states")
    pb$tick(0)
  }
   
  move <- function(state, depth = 0L) {
    if (depth > horizon)
      return()
    
    for (action in actions){
      next_states <- transition_function(model, action, state) 
        
      for (next_state in names(next_states)) {
        if (exists(next_state, envir = states))
          next()
          
        if (progress)
          pb$tick()
        
        assign(next_state, TRUE, envir = states)
        move(next_state, depth + 1L)
      }
    }
  }
  
  for (state in start_state) {
    assign(state, TRUE, envir = states)
    move(state)
  }
  
  ls(states)
}