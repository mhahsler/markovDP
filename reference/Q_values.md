# Q-Values

Several useful functions to deal with Q-values (action values) which map
each state/action pair to a utility value.

## Usage

``` r
Q_values(model, V = NULL, state = NULL)

Q_zero(model, value = 0)

Q_random(model, min = 1e-06, max = 1)
```

## Arguments

- model:

  an MDP problem specification.

- V:

  the state values. If `model` is a solved model, then the state values
  are taken from the solution.

- state:

  specify the state. If `NULL` then the Q-values for all states is
  returned as a matrix.

- value:

  value to initialize the Q-value matrix. Default is 0.

- min, max:

  range of the random values

## Value

`Q_values()` returns a state by action matrix specifying the Q-function,
i.e., the action value for executing each action in each state. The
Q-values are calculated from the value function (U) and the transition
model.

`Q_zero()` and `Q_random` return a matrix with q-values.

## Details

Implemented functions are:

- `Q_values()` gets the Q-values from the model. If the value function
  `V` is specified or if the policy of a solved model contains a value
  function, then the Q-values are approximated using the Bellman
  optimality equation:

  \$\$q\_\*(s,a) = \sum\_{s'} p(s'\|s,a) \[r(s,a,s') + \gamma
  v\_\*(s')\]\$\$

  Exact Q values are calculated if \\v = v\_\*\\, the optimal value
  function, otherwise we get an approximation that might not be
  consistent with \\v\\ or the implied policy. Q values can be used as
  the input for several other functions.

  For solvers, that directly calculate Q-values or an approximate
  Q-function, then the solver's values are returned.

- `Q_zero()` and `Q_random()` create initial Q value matrices for
  algorithms.

## References

Sutton, R. S., Barto, A. G. (2020). Reinforcement Learning: An
Introduction. Second edition. The MIT Press.

## See also

Other value_function:
[`bellman_update()`](http://michael.hahsler.net/markovDP/reference/bellman_update.md),
[`value_function()`](http://michael.hahsler.net/markovDP/reference/value_function.md)

## Author

Michael Hahsler

## Examples

``` r
data(Maze)
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

# create a random policy and calculate q-values
pi_random <- random_policy(Maze)
pi_random
#>     state  V action
#> 1  s(1,1) NA  right
#> 2  s(2,1) NA   left
#> 3  s(3,1) NA   down
#> 4  s(1,2) NA  right
#> 5  s(3,2) NA     up
#> 6  s(1,3) NA   down
#> 7  s(2,3) NA     up
#> 8  s(3,3) NA   left
#> 9  s(1,4) NA   left
#> 10 s(2,4) NA     up
#> 11 s(3,4) NA     up

V <- policy_evaluation(Maze, pi_random)
V
#>     s(1,1)     s(2,1)     s(3,1)     s(1,2)     s(3,2)     s(1,3)     s(2,3) 
#> -0.9509173 -4.8176833 -8.2898565 -0.4176471 -7.8997874 -0.3676471 -0.4823529 
#>     s(3,3)     s(1,4)     s(2,4)     s(3,4) 
#> -7.1191368  0.0000000  0.0000000 -1.7242440 

# calculate Q values
Q <- Q_values(Maze, V)
Q
#>                up      right       down       left
#> s(1,1) -0.9375903 -0.9509777 -4.0310031 -1.3775939
#> s(2,1) -1.7642705 -4.8182240 -7.6354218 -4.8182240
#> s(3,1) -5.5131111 -7.6705839 -8.2908495 -7.9826391
#> s(1,2) -0.5059741 -0.4176471 -0.5059741 -0.8842632
#> s(3,2) -7.9007292 -7.3152669 -7.9007292 -8.2518426
#> s(1,3) -0.2758824  0.6750000 -0.3676471 -0.4591176
#> s(2,3) -0.4823529 -1.5886784 -5.8835448 -1.1745607
#> s(3,3) -1.3882855 -2.1795442 -6.6977126 -7.1199789
#> s(1,4)  0.0000000  0.0000000  0.0000000  0.0000000
#> s(2,4)  0.0000000  0.0000000  0.0000000  0.0000000
#> s(3,4) -1.7243381 -1.6918196 -2.3037333 -6.0077339
```
