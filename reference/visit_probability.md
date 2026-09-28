# State Visit Probability

Calculates the state visit probability (the modified stationary
distribution) when following a policy from the start state.

## Usage

``` r
visit_probability(model, pi = NULL, start = NULL, method = "power", ...)
```

## Arguments

- model:

  a solved [MDP](http://michael.hahsler.net/markovDP/reference/MDP.md)
  object.

- pi:

  the used policy. If missing the policy in `model` is used.

- start:

  specification of the start distribution. If missing the specification
  in `model` is used.

- method:

  calculate the modified stationary distribution using `"power"` (power
  iteration) or `"sample"` (trajectory sampling).

- ...:

  further arguments are passed on to the internal implementations (see
  Details section).

## Value

a visit probability vector over all states.

## Details

The visit probability is the stationary distribution for the transition
matrix induced by the policy. To account for absorbing states, we modify
the transition matrix by setting all outgoing probabilities from
absorbing states to 0.

### Power iteration

The stationary distribution can be estimated as the sum of multiplying
the start distribution repeatedly with the modified transition matrix
induced by the policy. We stop multiplying when the largest difference
between entries in the two consecutive vectors is less then the extra
parameter:

- `min_err` stop criterion for the power iteration (default: `1e-6`).

- `sparse` logical; should a sparse transition matrix be used for the
  power iteration?

The resulting vector is normalized to probabilities.

### Sample method

The stationary distribution is calculated using `n` random walks. The
state visit counts are normalized to a probabilities. Additional
parameters are:

- `n` number of random walks (default `1000`).

- `horizon` maximal horizon used to stop a random walk if it has not
  reached an absorbing state.

## See also

Other policy:
[`action()`](http://michael.hahsler.net/markovDP/reference/action.md),
[`expected_return()`](http://michael.hahsler.net/markovDP/reference/expected_return.md),
[`greedy_action()`](http://michael.hahsler.net/markovDP/reference/greedy_action.md),
[`policy()`](http://michael.hahsler.net/markovDP/reference/policy.md),
[`policy_evaluation()`](http://michael.hahsler.net/markovDP/reference/policy_evaluation.md),
[`regret()`](http://michael.hahsler.net/markovDP/reference/regret.md)

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
visit_probability(sol)
#>      s(1,1)      s(2,1)      s(3,1)      s(1,2)      s(3,2)      s(1,3) 
#> 0.162710379 0.183049179 0.162710382 0.162710373 0.020338798 0.160481419 
#>      s(2,3)      s(3,3)      s(1,4)      s(2,4)      s(3,4) 
#> 0.017831260 0.000000000 0.128385085 0.001783124 0.000000000 

# gw_matrix also can calculate the visit_probability.
gw_matrix(sol, what = "visit_probability")
#>           [,1]      [,2]       [,3]        [,4]
#> [1,] 0.1627104 0.1627104 0.16048142 0.128385085
#> [2,] 0.1830492        NA 0.01783126 0.001783124
#> [3,] 0.1627104 0.0203388 0.00000000 0.000000000
```
