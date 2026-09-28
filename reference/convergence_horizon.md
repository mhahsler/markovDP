# Estimate the Convergence Horizon for an Infinite-Horizon MDP

Many sampling-based methods require a finite horizon. For infinite
horizons, discounting leads to convergences during a finite horizon.
This function estimates the number of steps till convergence using rules
of thumb.

## Usage

``` r
convergence_horizon(model, delta = 0.001, n_updates = 10)
```

## Arguments

- model:

  an MDP model.

- delta:

  maximum update error.

- n_updates:

  integer; average number of time each state is updated.

## Value

An estimated convergence horizon.

## Details

The horizon is estimated differently for the discounted and the
undiscounted case.

### Discounted Case

The effect of the largest reward \\R\_{\mathrm{max}}\\ update decreases
with \\t\\ as \\\delta_t = \gamma^t R\_{\mathrm{max}}\\. The convergence
horizon is estimated as the smallest \\t\\ for which \\\delta_t \<
\delta\\.

### Undiscounted Case

For the undiscounted case, episodes end when an absorbing state is
reached. It cannot be guaranteed that a model will reach an absorbing
state. To avoid infinite loops, we set the maximum horizon such that
each entry in the Q-table is on average updated `n_updates` times. This
is a very rough rule ot thumb.

## See also

Other solver:
[`schedule`](http://michael.hahsler.net/markovDP/reference/schedule.md),
[`solve_MDP()`](http://michael.hahsler.net/markovDP/reference/solve_MDP.md),
[`solve_MDP_APPROX()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_APPROX.md),
[`solve_MDP_DP()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_DP.md),
[`solve_MDP_LP()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_LP.md),
[`solve_MDP_MC()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_MC.md),
[`solve_MDP_PG()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_PG.md),
[`solve_MDP_SAMP()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_SAMP.md),
[`solve_MDP_TD()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_TD.md)

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

convergence_horizon(Maze)
#> Warning: discount needs to be <1 to guarantee convergence.
#>   Using a maximum horizon of |S| x |A| x n_updates = 440
#> [1] 440

# make the Maze into a discounted problem where future rewards count less.
Maze_discounted <- Maze
Maze_discounted$discount <- .9
Maze_discounted
#> MDPModel, MDP - Stuart Russell's 3x4 Maze
#>   Discount factor: 0.9
#>   Horizon: Inf epochs
#>   Size: 4 actions / 11 states
#>   Storage: transition prob as matrix / reward as matrix. Total size: 28.8 Kb
#>   Start: s(3,1)
#>   Model list components: ‘name’, ‘discount’, ‘horizon’, ‘states’,
#>     ‘actions’, ‘start’, ‘transition_model’, ‘reward’, ‘info’,
#>     ‘absorbing_states’

convergence_horizon(Maze_discounted)
#> [1] 66
```
