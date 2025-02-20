#' Find Reachable States 
#'
#' Finds the reachable state space from a MDP or MDPTF.
#'
#' There are three application cases for finding reachable states.
#' 
#' * For an MDP that has a state space defined, not all specified states might 
#'   be reachable. The function performs a (depth-limited) depth-first traversal 
#'   of the state space and returns a vector with the names of all 
#'   encountered states. This is used for example for [unreachable_states()]. 
#' 
#' * We may only have a transition model function following the specifications
#'   of an R transition function for an [MDP] which returns a probability 
#'   distribution over state transitions. To create a complete MDP, 
#'   we can use the found reachable states to create a complete MDP object. This
#'   search also used depth-first search of the state space.
#' 
#' * To use tabular methods for [MDPTF]s which specify a transition function that 
#'   returns the reward and the next state, we need to also specify the state space.
#'   Since an MDPTF can have a stochastic transition model, trajectory sampling is used.
#'   The horizon and the number of trajectories `n` has to be specified. **Note** that
#'   not all reachable states may be returned if some states have a very low probability 
#'   to be in the sample trajectories.
#'
#' @family MDP
#' @family MDPTF
#'
#' @author Michael Hahsler
#' 
#' @examples
#' # Example 1: Find the states of a simple MFPTD
#' 
#' line_maze <- MDPTF(actions = c("left", "right"), 
#'       start = s(0), 
#'       absorbing_states = rbind(s(-10), s(10)), 
#'       transition_func = function(model, state, action) {
#'             if (state == s(-10) || state == s(+10)) {
#'               return(list(reward = 0, state_prime = state))
#'             }
#'             
#'             
#'             reward <- 0
#'             if (action == "left") state <- state - 1
#'             if (action == "right") state <- state + 1
#'               
#'             if (state == s(-10) || state == s(+10)) {
#'                 reward <- 10
#'             }
#'         
#'             return(list(reward = reward, state_prime = state))
#'       },
#'       name = "line maze [-10,+10]"
#'       )
#'
#' # this model has no state specified
#' S(line_maze)
#' 
#' # find the states
#' states <- reachable_states(line_maze, horizon = 100)
#' states 
#' 
#' # set the states in the model
#' line_maze$states <- states
#' sol <- solve_MDP(line_maze, method = "TD:q_learning", 
#'   horizon = 100, n = 100, epsilon = .8)
#' 
#' policy(sol)
#' plot_value_function(sol)
#' 
#' # Example 2: Find the states to define a MDP for Tic-Tac-Toe
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
#'   
#'   state[action] <- player
#'   state
#' }
#' 
#' ttt_terminal <- function(state) {
#'   if (length(state) == 1L) state <- ttt_label2state(state)
#'   
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
#' # define the transition function: 
#' #     * return a probability vector for an action in a start state
#' #     * we define the special states 'win', 'loss', and 'draw'
#' P <- function(model, action, start.state) {
#'   action <- as.integer(action)
#'   
#'   # absorbing states
#'   if (start.state %in% c('win', 'loss', 'draw', 'illegal')) {
#'     return(structure(1, names = start.state))
#'   }
#'   
#'   # avoid illegal action by going to the very expensive illegal state
#'   if (!(action %in% ttt_available_actions(start.state))) {
#'     return(structure(1, names = "illegal"))
#'   }
#'   
#'   # make x's move
#'   next_state <- ttt_result(start.state, 'x', action)
#'   
#'   # terminal?
#'   term <- ttt_terminal(next_state)
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
#'   # fix terminal states
#'   term <- sapply(possible_end_states, ttt_terminal)
#'   possible_end_states <- sapply(possible_end_states, ttt_state2label)
#'   possible_end_states[term == 'x'] <- 'win'
#'   possible_end_states[term == 'o'] <- 'loss'
#'   possible_end_states[term == 'd'] <- 'draw'
#'   
#'   possible_end_states <- unique(possible_end_states)
#'   
#'   return(structure(rep(1 / length(possible_end_states), 
#'                       length(possible_end_states)), 
#'                    names = possible_end_states))
#' }
#' 
#' # define the reward
#' R <- rbind(
#'   R_(                    value = 0),
#'   R_(end.state = 'win',  value = +1),
#'   R_(end.state = 'loss', value = -1),
#'   R_(end.state = 'draw', value = +.5),
#'   R_(end.state = 'illegal', value = -Inf),
#'   # Note: there is no more reward once the agent is in a terminal state
#'   R_(start.state = 'win',  value = 0),
#'   R_(start.state = 'loss', value = 0),
#'   R_(start.state = 'draw', value = 0),
#'   R_(start.state = 'illegal', value = 0)
#' )
#' 
#' # start state
#' start <- ttt_state2label(ttt_empty_board())
#' start
#' 
#' # find the reachable state space
#' S <- union(c('win', 'loss', 'draw', 'illegal'),
#'           reachable_states(P, start_state = start, actions = A))
#' head(S)
#' 
#' tictactoe <- MDP(S, A, P, R, discount = 1, start = start, name = "TicTacToe")
#' tictactoe
#' 
#' # this MDP takes a about 30 seconds to solve using value iteration
#' # sol <- solve_MDP(tictactoe)
#' # policy(sol)[1:10, ]
#' 
#' @param model a MDP, MDPE or a MDP transition function.
#' @param progress logical; show a progress bar?
#' @param ... further arguments are passed on (e.g., to [`sample_MDP()`])
#' 
#' @returns a character vector with all reachable states.
#' @export
reachable_states <- function(model,
                             ...,
                             progress = TRUE) {
  UseMethod("reachable_states")
}


# use depth-first search
#' @importFrom fastmap fastmap faststack
#' @rdname reachable_states
#' @param horizon only return states reachable in the given horizon.
#' @export
reachable_states.MDP <- function(model,
                              horizon = Inf,
                              ...,
                              progress = TRUE) {
  reached <- fastmap()
  frontier <- faststack()
  
  for (start_state in start_vector(model, sparse = "states")) {
    frontier$push(start_state)
    # key: state label; value: depth
    reached$set(start_state, 0L)
  }
  
  if (progress) {
    pb <- my_progress_bar(N = length(S(model)), name = "unreachable_states")
    pb$tick(0)
  }
  
  while (frontier$size() > 0) {
    state <- frontier$pop()
    
    if (progress) {
      pb$tick()
    }
    
    # available_actions is slow!
    #for (action in available_actions(model, state)){
    for (action in A(model)) {
      next_states <- transition_matrix(model, action, state, sparse = "states")
      
      for (next_state in next_states) {
        if (reached$has(next_state))
          next()
        
        depth <- reached$get(state)
        if (depth >= horizon)
          next()
        
        frontier$push(next_state)
        reached$set(next_state, depth + 1L)
      }
    }
  }
  
  if (progress)
    pb$terminate()
  
  names(reached$as_list())
}

#' @rdname reachable_states
#' @param n number if sampled trajectories.
#' @export
reachable_states.MDPTF <- function(model,
                                   n = 100,
                                   horizon = NULL,
                                   ...,
                                   progress = TRUE) {
  horizon <- horizon %||% model$horizon
  
  samp <- sample_MDP(model, n = n, horizon = horizon, trajectories = FALSE, progress = progress, ...)
  
  sort(names(samp$state_cnt))
}


# model is a transition function for an MDP
#' @rdname reachable_states
#' @param actions labels of the available actions.
#' @param start_state label of the start state.
#' @export
reachable_states.function <- function(model,
                                   actions,
                                   start_state,
                                   horizon = Inf,
                                   ...,
                                   progress = TRUE) {
  transition_function <- model
    
  if (is.null(model)) {
      model <- list(actions = actions, 
                    start = start_state, 
                    transition_prob = transition_function)
    }
    
    if (progress) {
      pb <- my_progress_spinner(name = "find_reachable_states", ticks = "reached states")
      pb$tick(0)
    }
    
    # recursive depth-first search
    move <- function(state, depth = 0L) {
      if (depth > horizon)
        return()
      
      for (action in actions){
        next_states <- transition_function(model, action, state) 
        
        for (next_state in names(next_states)) {
          if (states$has(next_state))
            next()
          
          if (progress)
            pb$tick()
          
          states$set(next_state, TRUE)
          move(next_state, depth + 1L)
        }
      }
    }
    
    states <- fastmap()
    depth <- 0L
    
    
    for (state in start_state) {
      states$set(state, TRUE)
      move(state)
    }
    
    names(states$as_list())
  }
    
