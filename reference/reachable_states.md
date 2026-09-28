# Find Reachable States

Finds the reachable state space from an MDP or MDPSample.

## Usage

``` r
reachable_states(model, ..., progress = TRUE)

# S3 method for class 'MDPModel'
reachable_states(model, horizon = Inf, ..., progress = TRUE)

# S3 method for class 'MDPSample'
reachable_states(model, n = 100, horizon = NULL, ..., progress = TRUE)

# S3 method for class '`function`'
reachable_states(
  model,
  actions,
  start_state,
  horizon = Inf,
  ...,
  progress = TRUE
)
```

## Arguments

- model:

  an MDP or an MDP transition function.

- ...:

  further arguments are passed on (e.g., to
  [`sample_MDP()`](http://michael.hahsler.net/markovDP/reference/sample_MDP.md))

- progress:

  logical; show a progress bar?

- horizon:

  only return states reachable in the given horizon.

- n:

  number if sampled trajectories.

- actions:

  labels of the available actions.

- start_state:

  label of the start state.

## Value

a character vector with all reachable states.

## Details

There are three application cases for finding reachable states.

- For an MDP that has a state space defined, not all specified states
  might be reachable. The function performs a (depth-limited)
  depth-first traversal of the state space and returns a vector with the
  names of all encountered states. This is used for example for
  [`unreachable_states()`](http://michael.hahsler.net/markovDP/reference/unreachable_states.md).

- We may only have a transition model function following the
  specifications of an R transition function for an
  [MDPModel](http://michael.hahsler.net/markovDP/reference/MDP.md) which
  returns a probability distribution over state transitions. To create a
  complete MDP, we can use the found reachable states to create a
  complete MDP object. This search also used depth-first search of the
  state space.

- To use tabular methods for
  [MDPSample](http://michael.hahsler.net/markovDP/reference/MDPSample.md)
  which specify a transition function that returns the reward and the
  next state, we need to also specify the state space. Since an
  MDPSample can have a stochastic transition model, trajectory sampling
  is used. The horizon and the number of trajectories `n` has to be
  specified.

  **Notes:**

  - Not all reachable states may be returned if some states have a very
    low probability to be in the sample trajectories.

  - A finite horizon is needed for sampling. By default the model
    horizon is used. Infinite horizon is capped at 1000. Larger values
    can be specified manually.

## See also

Other MDP:
[`MDP()`](http://michael.hahsler.net/markovDP/reference/MDP.md),
[`absorbing_states()`](http://michael.hahsler.net/markovDP/reference/absorbing_states.md),
[`act()`](http://michael.hahsler.net/markovDP/reference/act.md),
[`action_state_helpers`](http://michael.hahsler.net/markovDP/reference/action_state_helpers.md),
[`available_actions()`](http://michael.hahsler.net/markovDP/reference/available_actions.md),
[`find_reachable_states()`](http://michael.hahsler.net/markovDP/reference/find_reachable_states.md),
[`sample_MDP()`](http://michael.hahsler.net/markovDP/reference/sample_MDP.md),
[`sample_MDP.MDPSample()`](http://michael.hahsler.net/markovDP/reference/sample_MDP.MDPSample.md),
[`start`](http://michael.hahsler.net/markovDP/reference/start.md),
[`transition_graph()`](http://michael.hahsler.net/markovDP/reference/transition_graph.md),
[`transition_matrix()`](http://michael.hahsler.net/markovDP/reference/accessors.md),
[`unreachable_states()`](http://michael.hahsler.net/markovDP/reference/unreachable_states.md)

Other MDPSample:
[`MDPSample()`](http://michael.hahsler.net/markovDP/reference/MDPSample.md),
[`absorbing_states()`](http://michael.hahsler.net/markovDP/reference/absorbing_states.md),
[`act()`](http://michael.hahsler.net/markovDP/reference/act.md),
[`action_state_helpers`](http://michael.hahsler.net/markovDP/reference/action_state_helpers.md),
[`sample_MDP.MDPSample()`](http://michael.hahsler.net/markovDP/reference/sample_MDP.MDPSample.md),
[`solve_MDP_PG()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_PG.md),
[`start`](http://michael.hahsler.net/markovDP/reference/start.md)

## Author

Michael Hahsler

## Examples

``` r
# Example 1: Find the reachable states of a simple MDP with only sample access

line_maze <- MDPSample(actions = c("left", "right"), 
      start = s(0), 
      absorbing_states = rbind(s(-10), s(10)), 
      transition_model = function(model, action, state) {
            if (state == s(-10) || state == s(+10)) {
              return(list(reward = 0, state_prime = state))
            }
            
            reward <- 0
            if (action == "left") state <- state - 1
            if (action == "right") state <- state + 1
              
            if (state == s(-10) || state == s(+10)) {
                reward <- 10
            }
        
            return(list(reward = reward, state_prime = state))
      },
      name = "line maze [-10,+10]"
      )

# this model has no state specified
line_maze
#> MDPSample, MDP - line maze [-10,+10]
#>   Discount factor: 0.9
#>   Horizon: Inf epochs
#>   Size: 2 actions / undefined state space
#>   Start: s(0)
#>   List components: ‘name’, ‘discount’, ‘horizon’, ‘actions’, ‘states’,
#>     ‘start’, ‘absorbing_states’, ‘transition_model’, ‘info’
S(line_maze)
#> NULL

# find the states
states <- reachable_states(line_maze, horizon = 100)
states 
#>  [1] "s(-1)"  "s(-10)" "s(-2)"  "s(-3)"  "s(-4)"  "s(-5)"  "s(-6)"  "s(-7)" 
#>  [9] "s(-8)"  "s(-9)"  "s(0)"   "s(1)"   "s(10)"  "s(2)"   "s(3)"   "s(4)"  
#> [17] "s(5)"   "s(6)"   "s(7)"   "s(8)"   "s(9)"  

# set the states in the model
line_maze$states <- states
line_maze
#> MDPSample, MDP - line maze [-10,+10]
#>   Discount factor: 0.9
#>   Horizon: Inf epochs
#>   Size: 2 actions / 21 states
#>   Start: s(0)
#>   List components: ‘name’, ‘discount’, ‘horizon’, ‘actions’, ‘states’,
#>     ‘start’, ‘absorbing_states’, ‘transition_model’, ‘info’

sol <- solve_MDP(line_maze, method = "TD:q_learning", 
  horizon = 100, n = 100, epsilon = .8)

policy(sol)
#>     state         V action
#> 1   s(-1) 2.5001888   left
#> 2  s(-10) 0.0000000  right
#> 3   s(-2) 3.2589494   left
#> 4   s(-3) 4.1230376   left
#> 5   s(-4) 5.0752757   left
#> 6   s(-5) 6.0703161   left
#> 7   s(-6) 7.0287817   left
#> 8   s(-7) 7.9976524   left
#> 9   s(-8) 8.9895388   left
#> 10  s(-9) 9.9990328   left
#> 11   s(0) 1.7037361   left
#> 12   s(1) 1.3041920   left
#> 13  s(10) 0.0000000  right
#> 14   s(2) 1.0743627   left
#> 15   s(3) 0.8032204   left
#> 16   s(4) 0.3581580   left
#> 17   s(5) 0.1192231   left
#> 18   s(6) 0.1333712  right
#> 19   s(7) 0.5072045  right
#> 20   s(8) 1.8258678  right
#> 21   s(9) 5.7824847  right
plot_value_function(sol)


# Example 2: Find the states to define an MDP for Tic-Tac-Toe

# state description: matrix with the characters _, x, and o
#                    can be converted into a label of 9 characters

# set of actions
A <- as.character(1:9)

# helper functions
ttt_empty_board <- function() matrix('_', ncol = 3, nrow = 3)

ttt_state2label <- function(state) paste(state, collapse = '')

ttt_label2state <- function(label) matrix(strsplit(label, "")[[1]], 
                                          nrow = 3, ncol = 3)

ttt_available_actions <- function(state) {
  if (length(state) == 1L) state <- ttt_label2state(state)
  which(state == "_")
}

ttt_result <- function(state, player, action) {
  if (length(state) == 1L) state <- ttt_label2state(state)
  
  if (state[action] != "_")
    stop("Illegal action.")
  
  state[action] <- player
  state
}

ttt_terminal <- function(state) {
  if (length(state) == 1L) state <- ttt_label2state(state)
  
  # Check the board for a win and return one of 
  # 'x', 'o', 'd' (draw), or 'n' (for next move)
  win_possibilities <- rbind(state, 
                             t(state), 
                             diag(state), 
                             diag(t(state)))

  wins <- apply(win_possibilities, MARGIN = 1, FUN = function(x) {
    if (x[1] != '_' && length(unique(x)) == 1) x[1]
    else '_'
  })

  if (any(wins == 'x')) 
    return('x')

  if (any(wins == 'o')) 
    return('o')

  # Check for draw
  if (sum(state == '_') < 1)
    return('d')

  return('n')
}

# define the transition function: 
#     * return a probability vector for an action in a start state
#     * we define the special states 'win', 'loss', and 'draw'
P <- function(model, action, start.state) {
  action <- as.integer(action)
  
  # absorbing states
  if (start.state %in% c('win', 'loss', 'draw', 'illegal')) {
    return(structure(1, names = start.state))
  }
  
  # avoid illegal action by going to the very expensive illegal state
  if (!(action %in% ttt_available_actions(start.state))) {
    return(structure(1, names = "illegal"))
  }
  
  # make x's move
  next_state <- ttt_result(start.state, 'x', action)
  
  # terminal?
  term <- ttt_terminal(next_state)
  if (term == 'x') {
    return(structure(1, names = "win"))
  } else if (term == 'o') {
    return(structure(1, names = "loss"))
  } else if (term == 'd') {
    return(structure(1, names = "draw"))
  }
  
  # it is o's turn
  actions_of_o <- ttt_available_actions(next_state)
  possible_end_states <- lapply(
    actions_of_o,
    FUN = function(a)
      ttt_result(next_state, 'o', a)
  )
  
  # fix terminal states
  term <- sapply(possible_end_states, ttt_terminal)
  possible_end_states <- sapply(possible_end_states, ttt_state2label)
  possible_end_states[term == 'x'] <- 'win'
  possible_end_states[term == 'o'] <- 'loss'
  possible_end_states[term == 'd'] <- 'draw'
  
  possible_end_states <- unique(possible_end_states)
  
  return(structure(rep(1 / length(possible_end_states), 
                      length(possible_end_states)), 
                   names = possible_end_states))
}

# define the reward
R <- rbind(
  R_(                    value = 0),
  R_(end.state = 'win',  value = +1),
  R_(end.state = 'loss', value = -1),
  R_(end.state = 'draw', value = +.5),
  R_(end.state = 'illegal', value = -Inf),
  # Note: there is no more reward once the agent is in a terminal state
  R_(start.state = 'win',  value = 0),
  R_(start.state = 'loss', value = 0),
  R_(start.state = 'draw', value = 0),
  R_(start.state = 'illegal', value = 0)
)

# start state
start <- ttt_state2label(ttt_empty_board())
start
#> [1] "_________"

# find the reachable state space
S <- union(c('win', 'loss', 'draw', 'illegal'),
          reachable_states(P, start_state = start, actions = A))
head(S)
#> [1] "win"       "loss"      "draw"      "illegal"   "xoxo_ox__" "xoxoo_xxo"

tictactoe <- MDP(S, A, P, R, discount = 1, start = start, name = "TicTacToe")
tictactoe
#> MDPModel, MDP - TicTacToe
#>   Discount factor: 1
#>   Horizon: Inf epochs
#>   Size: 9 actions / 2527 states
#>   Storage: transition prob as function / reward as data.frame. Total size: 201.3 Kb
#>   Start: _________
#>   Model list components: ‘name’, ‘discount’, ‘horizon’, ‘states’,
#>     ‘actions’, ‘start’, ‘transition_model’, ‘reward’, ‘info’

# this MDP takes about 30 seconds to solve using value iteration
# sol <- solve_MDP(tictactoe)
# policy(sol)[1:10, ]
```
