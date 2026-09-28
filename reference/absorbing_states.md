# Absorbing States

Find absorbing states using the transition model.

## Usage

``` r
absorbing_states(model, state = NULL, ...)

# S3 method for class 'MDP'
absorbing_states(
  model,
  state = NULL,
  sparse = "states",
  use_precomputed = TRUE,
  ...
)

# S3 method for class 'MDPSample'
absorbing_states(
  model,
  state = NULL,
  sparse = "features",
  use_precomputed = TRUE,
  ...
)
```

## Arguments

- model:

  a [MDP](http://michael.hahsler.net/markovDP/reference/MDP.md) object.

- state:

  a single state to check. This can be much faster if the model contains
  a transition model implemented as a function. `NULL` means all states
  are checked.

- ...:

  further arguments are passed on.

- sparse:

  logical; if return a sparse logical vector?

- use_precomputed:

  logical; should precomputed values in the MDP be used?

## Value

`absorbing_states()` returns a logical vector indicating if the states
are absorbing (terminal).

## Details

The function `absorbing_states()` checks if a state or a set of states
are absorbing (terminal states). A state is absorbing if there is for
all actions a probability of 1 for staying in the state.

## See also

Other MDP:
[`MDP()`](http://michael.hahsler.net/markovDP/reference/MDP.md),
[`act()`](http://michael.hahsler.net/markovDP/reference/act.md),
[`action_state_helpers`](http://michael.hahsler.net/markovDP/reference/action_state_helpers.md),
[`available_actions()`](http://michael.hahsler.net/markovDP/reference/available_actions.md),
[`find_reachable_states()`](http://michael.hahsler.net/markovDP/reference/find_reachable_states.md),
[`reachable_states()`](http://michael.hahsler.net/markovDP/reference/reachable_states.md),
[`sample_MDP()`](http://michael.hahsler.net/markovDP/reference/sample_MDP.md),
[`sample_MDP.MDPSample()`](http://michael.hahsler.net/markovDP/reference/sample_MDP.MDPSample.md),
[`start`](http://michael.hahsler.net/markovDP/reference/start.md),
[`transition_graph()`](http://michael.hahsler.net/markovDP/reference/transition_graph.md),
[`transition_matrix()`](http://michael.hahsler.net/markovDP/reference/accessors.md),
[`unreachable_states()`](http://michael.hahsler.net/markovDP/reference/unreachable_states.md)

Other MDPSample:
[`MDPSample()`](http://michael.hahsler.net/markovDP/reference/MDPSample.md),
[`act()`](http://michael.hahsler.net/markovDP/reference/act.md),
[`action_state_helpers`](http://michael.hahsler.net/markovDP/reference/action_state_helpers.md),
[`reachable_states()`](http://michael.hahsler.net/markovDP/reference/reachable_states.md),
[`sample_MDP.MDPSample()`](http://michael.hahsler.net/markovDP/reference/sample_MDP.MDPSample.md),
[`solve_MDP_PG()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_PG.md),
[`start`](http://michael.hahsler.net/markovDP/reference/start.md)

## Author

Michael Hahsler

## Examples

``` r
data(Maze)

gw_matrix(Maze)
#>      [,1]     [,2]     [,3]     [,4]    
#> [1,] "s(1,1)" "s(1,2)" "s(1,3)" "s(1,4)"
#> [2,] "s(2,1)" NA       "s(2,3)" "s(2,4)"
#> [3,] "s(3,1)" "s(3,2)" "s(3,3)" "s(3,4)"
gw_matrix(Maze, what = "labels")
#>      [,1]    [,2] [,3] [,4]      
#> [1,] ""      ""   ""   "Goal: +1"
#> [2,] ""      "X"  ""   "-1"      
#> [3,] "Start" ""   ""   ""        
gw_matrix(Maze, what = "absorbing")
#>       [,1]  [,2]  [,3]  [,4]
#> [1,] FALSE FALSE FALSE  TRUE
#> [2,] FALSE    NA FALSE  TRUE
#> [3,] FALSE FALSE FALSE FALSE

# -1 and +1 are absorbing states
absorbing_states(Maze)
#> [1] "s(1,4)" "s(2,4)"
absorbing_states(Maze, sparse = FALSE)
#> s(1,1) s(2,1) s(3,1) s(1,2) s(3,2) s(1,3) s(2,3) s(3,3) s(1,4) s(2,4) s(3,4) 
#>  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE   TRUE   TRUE  FALSE 
absorbing_states(Maze, sparse = "states")
#> [1] "s(1,4)" "s(2,4)"

# check individual states
absorbing_states(Maze, "s(1,1)")
#> [1] FALSE
absorbing_states(Maze, "s(1,4)")
#> [1] TRUE
```
