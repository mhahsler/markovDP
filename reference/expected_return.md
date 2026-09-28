# Calculate the Expected Return of a Policy

This function calculates the expected total return for an MDP policy
given a start state (distribution). The value is calculated using the
value function stored in the MDP solution.

## Usage

``` r
expected_return(model, ...)

# S3 method for class 'MDP'
expected_return(model, start = NULL, method = "solution", ...)
```

## Arguments

- model:

  a solved [MDP](http://michael.hahsler.net/markovDP/reference/MDP.md)
  object.

- ...:

  further arguments are passed on to
  [`policy_evaluation()`](http://michael.hahsler.net/markovDP/reference/policy_evaluation.md)
  or
  [`sample_MDP()`](http://michael.hahsler.net/markovDP/reference/sample_MDP.md).

- start:

  specification of the current state (see argument start in
  [MDP](http://michael.hahsler.net/markovDP/reference/MDP.md) for
  details). By default the start state defined in the model as start is
  used. Multiple states can be specified as rows in a matrix.

- method:

  `"solution"` uses the converged value function stored in the solved
  model, `"policy_evaluation"` estimates the value function, and
  `"sample"` calculates the average return by sampling episodes from the
  model.

## Value

`expected_return()` returns a vector of returns, one for each start
state if a matrix is specified.

- state:

  start state to calculate the return for. If `NULL` then the start
  state of model is used.

## Details

The return is typically calculated using the value function of the
solution. If these are not available, then
[`sample_MDP()`](http://michael.hahsler.net/markovDP/reference/sample_MDP.md)
is used instead with a warning.

## See also

Other policy:
[`action()`](http://michael.hahsler.net/markovDP/reference/action.md),
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
gw_matrix(Maze)
#>      [,1]     [,2]     [,3]     [,4]    
#> [1,] "s(1,1)" "s(1,2)" "s(1,3)" "s(1,4)"
#> [2,] "s(2,1)" NA       "s(2,3)" "s(2,4)"
#> [3,] "s(3,1)" "s(3,2)" "s(3,3)" "s(3,4)"

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
#> 10 s(2,4) 0.0000000     up
#> 11 s(3,4) 0.3872455   left

# return for the start state s(3,1) specified in the model
expected_return(sol)
#> [1] 0.7052527

# return for starting next to the goal at s(1,3)
expected_return(sol, start = "s(1,3)")
#> [1] 0.9178082

# expected return when we start from a random state as returned from the solver
expected_return(sol, start = "uniform")
#> [1] 0.5797937

# estimate the return using sampling following the policy
expected_return(sol, method = "sample", start = "uniform", n = 10000, horizon = 1000)
#> [1] 0.574804
```
