# Perform an Action

Performs an action in a state and returns the new state and reward.

## Usage

``` r
act(model, state, action, fast = FALSE, ...)

# S3 method for class 'MDPModel'
act(model, state, action = NULL, fast = FALSE, ...)

# S3 method for class 'MDPSample'
act(model, state, action, fast = FALSE, ...)
```

## Arguments

- model:

  an MDP model.

- state:

  the current state.

- action:

  the chosen action. If the action is not specified (`NULL`) and the MDP
  model contains a policy, then the action is chosen according to the
  policy.

- fast:

  logical; if `TRUE` then extra state id to label conversions are
  avoided.

- ...:

  if action is unspecified, then the additional parameters are passed on
  to
  [`action()`](http://michael.hahsler.net/markovDP/reference/action.md)
  to determine the action using the model's policy.

## Value

a named list with the `reward` and the next state `state_prime`.

## See also

Other MDP:
[`MDP()`](http://michael.hahsler.net/markovDP/reference/MDP.md),
[`absorbing_states()`](http://michael.hahsler.net/markovDP/reference/absorbing_states.md),
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
[`absorbing_states()`](http://michael.hahsler.net/markovDP/reference/absorbing_states.md),
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

act(Maze, "s(1,3)", "right")
#> $reward
#> [1] 0.96
#> 
#> $state_prime
#> [1] s(1,4)
#> 11 Levels: s(1,1) s(2,1) s(3,1) s(1,2) s(3,2) s(1,3) s(2,3) s(3,3) ... s(3,4)
#> 

# solve the maze and then ask for actions using the policy
sol <- solve_MDP(Maze)
act(sol, "s(1,3)")
#> $reward
#> [1] 0.96
#> 
#> $state_prime
#> [1] s(1,4)
#> 11 Levels: s(1,1) s(2,1) s(3,1) s(1,2) s(3,2) s(1,3) s(2,3) s(3,3) ... s(3,4)
#> 

# make the policy in sol epsilon-soft and ask 10 times for the action
replicate(10, act(sol, "s(1,3)", epsilon = .2))
#>             [,1]   [,2]   [,3]   [,4]   [,5]   [,6]   [,7]   [,8]   [,9]  
#> reward      0.96   -0.04  0.96   0.96   0.96   0.96   -0.04  -0.04  0.96  
#> state_prime s(1,4) s(2,3) s(1,4) s(1,4) s(1,4) s(1,4) s(2,3) s(1,3) s(1,4)
#>             [,10] 
#> reward      -0.04 
#> state_prime s(1,3)
```
