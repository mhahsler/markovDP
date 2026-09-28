# Choose an Action Given a Policy

Returns an action given a deterministic policy. The policy can be made
epsilon-soft.

## Usage

``` r
action(model, state, epsilon = 0, epoch = 1, as = "factor", ...)
```

## Arguments

- model:

  a solved [MDP](http://michael.hahsler.net/markovDP/reference/MDP.md).

- state:

  the state.

- epsilon:

  make the policy epsilon soft.

- epoch:

  what epoch of the policy should be used. Use 1 for converged policies.

- as:

  string, format for returning the action (e.g., `"factor"`, `"id"`,
  `"label"`).

- ...:

  further parameters are passed on.

## Value

The name of the optimal action as a factor.

## See also

Other policy:
[`expected_return()`](http://michael.hahsler.net/markovDP/reference/expected_return.md),
[`greedy_action()`](http://michael.hahsler.net/markovDP/reference/greedy_action.md),
[`policy()`](http://michael.hahsler.net/markovDP/reference/policy.md),
[`policy_evaluation()`](http://michael.hahsler.net/markovDP/reference/policy_evaluation.md),
[`regret()`](http://michael.hahsler.net/markovDP/reference/regret.md),
[`visit_probability()`](http://michael.hahsler.net/markovDP/reference/visit_probability.md)

## Author

Michael Hahsler

## Examples

``` r
data("Maze")
Maze
#> MDPModel, MDP - Stuart Russell's 3x4 Maze
#>   Discount factor: 1
#>   Horizon: Inf epochs
#>   Size: 4 actions / 11 states
#>   Storage: transition prob as matrix / reward as matrix. Total size: 28.8 Kb
#>   Start: s(3,1)
#>   Model list components: ‘name’, ‘discount’, ‘horizon’, ‘states’,
#>     ‘actions’, ‘start’, ‘transition_model’, ‘reward’, ‘info’,
#>     ‘absorbing_states’

sol <- solve_MDP(Maze)
policy(sol)
#>     state         V action
#> 1  s(1,1) 0.8115564  right
#> 2  s(2,1) 0.7615521     up
#> 3  s(3,1) 0.7052527     up
#> 4  s(1,2) 0.8678082  right
#> 5  s(3,2) 0.6551492   left
#> 6  s(1,3) 0.9178082  right
#> 7  s(2,3) 0.6602740     up
#> 8  s(3,3) 0.6110843   left
#> 9  s(1,4) 0.0000000   down
#> 10 s(2,4) 0.0000000   left
#> 11 s(3,4) 0.3872455   left

action(sol, state = "s(1,3)")
#> [1] right
#> Levels: up right down left

## choose from an epsilon-soft policy
table(replicate(100, action(sol, state = "s(1,3)", epsilon = 0.1)))
#> 
#>    up right  down  left 
#>     1    96     2     1 
```
