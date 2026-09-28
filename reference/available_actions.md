# Available Actions in a State

Determine the set of actions available in a state.

## Usage

``` r
available_actions(model, state, neg_inf_reward = TRUE, stay_in_place = FALSE)
```

## Arguments

- model:

  a [MDP](http://michael.hahsler.net/markovDP/reference/MDP.md) object.

- state:

  a character vector specifying the states.

- neg_inf_reward:

  logical; consider an action that produced `-Inf` reward to all end
  states unavailable?

- stay_in_place:

  logical; consider an action that results in the same state with a
  probability of 1 as unavailable. Note that this means that absorbing
  states have no available action!

## Value

a character vector with the available actions.

a vector with the available actions.

## Details

Unavailable actions are modeled as actions that have an immediate reward
of `-Inf` in the reward function. For a maze, also actions that do not
change the state can be considered unavailable.

## See also

Other MDP:
[`MDP()`](http://michael.hahsler.net/markovDP/reference/MDP.md),
[`absorbing_states()`](http://michael.hahsler.net/markovDP/reference/absorbing_states.md),
[`act()`](http://michael.hahsler.net/markovDP/reference/act.md),
[`action_state_helpers`](http://michael.hahsler.net/markovDP/reference/action_state_helpers.md),
[`find_reachable_states()`](http://michael.hahsler.net/markovDP/reference/find_reachable_states.md),
[`reachable_states()`](http://michael.hahsler.net/markovDP/reference/reachable_states.md),
[`sample_MDP()`](http://michael.hahsler.net/markovDP/reference/sample_MDP.md),
[`sample_MDP.MDPSample()`](http://michael.hahsler.net/markovDP/reference/sample_MDP.MDPSample.md),
[`start`](http://michael.hahsler.net/markovDP/reference/start.md),
[`transition_graph()`](http://michael.hahsler.net/markovDP/reference/transition_graph.md),
[`transition_matrix()`](http://michael.hahsler.net/markovDP/reference/accessors.md),
[`unreachable_states()`](http://michael.hahsler.net/markovDP/reference/unreachable_states.md)

## Author

Michael Hahsler

## Examples

``` r
data(DynaMaze)
gw_plot(DynaMaze)


# The following actions are always available:
DynaMaze$actions
#> [1] "up"    "right" "down"  "left" 

# only right and down is unavailable for s(1,1) because they
#   make the agent stay in place.
available_actions(DynaMaze, state = "s(1,1)", stay_in_place = TRUE)
#>           up right down  left
#> s(1,1) FALSE  TRUE TRUE FALSE

# An action that leaves the grid currently is allowed but does not do
# anything.
act(DynaMaze, "s(1,1)", "up")
#> $reward
#> [1] -1
#> 
#> $state_prime
#> [1] s(1,1)
#> 47 Levels: s(1,1) s(2,1) s(3,1) s(4,1) s(5,1) s(6,1) s(1,2) s(2,2) ... s(6,9)
#> 
```
