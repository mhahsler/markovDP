# Solve MDPs using Linear Programming

Solve discounted, infinite horizon MDPs via linear programming.

## Usage

``` r
solve_MDP_LP(
  model,
  method = "LP",
  horizon = NULL,
  discount = NULL,
  inf = 1000,
  lpSolve_args = list(),
  ...,
  matrix = NULL,
  continue = FALSE,
  verbose = FALSE,
  progress = NULL
)
```

## Arguments

- model:

  an MDP problem specification.

- method:

  string; one of the following solution methods: `'LP'`

- horizon:

  Only infinite-horizon MDPs with `horizon = Inf` are supported.

- discount:

  only undiscounted MDPs with `discount = 1` are supported.

- inf:

  value used for infinity when calling
  [`lpSolve::lp()`](https://rdrr.io/pkg/lpSolve/man/lp.html). This
  should me much larger than the largest absolute reward in the model.

- lpSolve_args:

  a list with additional arguments passed on to
  [`lpSolve::lp()`](https://rdrr.io/pkg/lpSolve/man/lp.html).

- ...:

  further parameters are passed on to the solver function.

- matrix:

  logical; if `TRUE` then matrices for the transition model and the
  reward function are taken from the model first. This can be slow if
  functions need to be converted or do not fit into memory if the models
  are large. If these components are already matrices, then this is very
  fast. For `FALSE`, the transition probabilities and the reward is
  extracted when needed. This is slower, but removes the time and memory
  requirements needed to calculate the matrices.

- continue:

  logical; show a progress bar with estimated time for completion.

- verbose:

  logical or a numeric verbose level; if set to `TRUE` or `1`, the
  function displays the used algorithm parameters and progress
  information. Levels `>1` provide more detailed solver output in the R
  console.

- progress:

  not supported by this solver.

## Value

[`solve_MDP()`](http://michael.hahsler.net/markovDP/reference/solve_MDP.md)
returns an object of class MDP or MDPSample which is a list with the
model specifications (`model`), the solution (`solution`). The solution
is a list with the elements that depend on the used method. Common
elements are:

- `method` with the name of the used method

- parameters used.

- `converged` did the algorithm converge (`NA`) for finite-horizon
  problems.

- `policy` a list representing the policy graph. The list only has one
  element for converged solutions.

## Details

A linear programming formulation was developed by Manne (1960) and
further described by Puterman (1996). For the optimal value function,
the Bellman equation holds:

\$\$ v^\*(s) = \max\_{a \in \mathcal{A}}\sum\_{s' \in \mathcal{S}} p(s,
a, s') \[ r(s, a, s') + \gamma v^\*(s')\]\\ \forall a\in \mathcal{A}, s
\in \mathcal{S} \$\$

The maximization problem can reformulate as a minimization with a linear
constraint for each state action pair. The optimal value function can be
found by solving the following linear program: \$\$\text{min}
\sum\_{s\in S} v(s)\$\$ subject to \$\$v(s) \ge \sum\_{s' \in
\mathcal{S}} p(s, a, s')\[r(s, a, s') + \gamma v(s')\],\\ \forall a\in
\mathcal{A}, s \in \mathcal{S} \$\$

Note:

- The discounting factor has to be strictly less than 1.

- The used solver does not support infinity and a sufficiently large
  value needs to be used instead (see parameter `inf`).

- Additional parameters to to `solve_MDP` are passed on to
  [`lpSolve::lp()`](https://rdrr.io/pkg/lpSolve/man/lp.html).

## References

Manne, Alan. 1960. "On the Job-Shop Scheduling Problem." Operations
Research 8 (2): 219-23.
[doi:10.1287/opre.8.2.219](https://doi.org/10.1287/opre.8.2.219) .

Puterman, Martin L. 1996. Markov decision processes: discrete stochastic
dynamic programming. John Wiley & Sons.

## See also

Other solver:
[`convergence_horizon()`](http://michael.hahsler.net/markovDP/reference/convergence_horizon.md),
[`schedule`](http://michael.hahsler.net/markovDP/reference/schedule.md),
[`solve_MDP()`](http://michael.hahsler.net/markovDP/reference/solve_MDP.md),
[`solve_MDP_APPROX()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_APPROX.md),
[`solve_MDP_DP()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_DP.md),
[`solve_MDP_MC()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_MC.md),
[`solve_MDP_PG()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_PG.md),
[`solve_MDP_SAMP()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_SAMP.md),
[`solve_MDP_TD()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_TD.md)

## Examples

``` r
data(Maze)

# we change the discount to 0.9 since LP is only implemented for discounted MDPs.  
maze_solved <- solve_MDP(Maze, discount = 0.9, method = "LP:LP", verbose = TRUE)
#> creating 44 constraints ... took 0.001 seconds.
#> running LP solver ... took 0.001 seconds.
maze_solved
#> MDPModel, MDP - Stuart Russell's 3x4 Maze
#>   Discount factor: 0.9
#>   Horizon: Inf epochs
#>   Size: 4 actions / 11 states
#>   Storage: transition prob as matrix / reward as matrix. Total size: 47.3 Kb
#>   Start: s(3,1)
#>   Model list components: ‘name’, ‘discount’, ‘horizon’, ‘states’,
#>     ‘actions’, ‘start’, ‘transition_model’, ‘reward’, ‘info’,
#>     ‘absorbing_states’, ‘solution’
#> 
#>   Solved:
#>     Method: ‘lp’
#>     Solution converged: TRUE
#>   Solution list components: ‘method’, ‘policy’, ‘converged’,
#>     ‘solver_out’
policy(maze_solved)
#>     state         V action
#> 1  s(1,1) 0.5810788  right
#> 2  s(2,1) 0.4614351     up
#> 3  s(3,1) 0.3508265     up
#> 4  s(1,2) 0.7322953  right
#> 5  s(3,2) 0.3002100  right
#> 6  s(1,3) 0.8895585  right
#> 7  s(2,3) 0.5499803     up
#> 8  s(3,3) 0.3974613     up
#> 9  s(1,4) 0.0000000   left
#> 10 s(2,4) 0.0000000     up
#> 11 s(3,4) 0.1606287   left
```
