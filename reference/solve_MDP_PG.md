# Solve MDPs with Policy Gradient Methods

Not implemented yet: Solve the MDP control problem using a parameterized
policy and policy gradient methods. The implemented method one-step
actor-critic control and actor-critic control with eligibility traces.

## Usage

``` r
solve_MDP_PG(
  model,
  method = "actor-critic",
  horizon = NULL,
  discount = NULL,
  alpha_actor = schedule_exp(0.2, 0.1),
  alpha_critic = schedule_exp(0.2, 0.1),
  epsilon = schedule_exp(1, 0.1),
  lambda = 0,
  n,
  w = NULL,
  theta = NULL,
  ...,
  matrix = TRUE,
  continue = FALSE,
  progress = TRUE,
  verbose = FALSE
)
```

## Arguments

- model:

  an MDP problem specification.

- method:

  string; one of the following solution methods: `'sarsa'`

- horizon:

  an integer with the number of epochs for problems with a finite
  planning horizon. If set to `Inf`, the algorithm continues running
  iterations till it converges to the infinite horizon solution. If
  `NULL`, then the horizon specified in `model` will be used.

- discount:

  discount factor in range \\(0, 1\]\\. If `NULL`, then the discount
  factor specified in `model` will be used.

- alpha_actor, alpha_critic:

  alpha schedules

- epsilon:

  used for the \\\epsilon\\-greedy behavior policies. A scalar value
  between 0 and 1 or a
  [schedule](http://michael.hahsler.net/markovDP/reference/schedule.md).

- lambda:

  the trace-decay parameter for the an accumulating trace. If
  `lambda = 0` then 1-step Sarsa is used.

- n:

  number of episodes used for learning.

- w:

  an initial weight vector. By default a vector with 0s is used.

- theta:

  parameter...

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

- progress:

  logical; show a progress bar with estimated time for completion.

- verbose:

  logical or a numeric verbose level; if set to `TRUE` or `1`, the
  function displays the used algorithm parameters and progress
  information. Levels `>1` provide more detailed solver output in the R
  console.

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

## References

Sutton, Richard S., and Andrew G. Barto. 2018. Reinforcement Learning:
An Introduction. Second. The MIT Press.
[http://incompleteideas.net/book/the-book-2nd.html](http://incompleteideas.net/book/the-book-2nd.md).

## See also

Other solver:
[`convergence_horizon()`](http://michael.hahsler.net/markovDP/reference/convergence_horizon.md),
[`schedule`](http://michael.hahsler.net/markovDP/reference/schedule.md),
[`solve_MDP()`](http://michael.hahsler.net/markovDP/reference/solve_MDP.md),
[`solve_MDP_APPROX()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_APPROX.md),
[`solve_MDP_DP()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_DP.md),
[`solve_MDP_LP()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_LP.md),
[`solve_MDP_MC()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_MC.md),
[`solve_MDP_SAMP()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_SAMP.md),
[`solve_MDP_TD()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_TD.md)

Other MDPSample:
[`MDPSample()`](http://michael.hahsler.net/markovDP/reference/MDPSample.md),
[`absorbing_states()`](http://michael.hahsler.net/markovDP/reference/absorbing_states.md),
[`act()`](http://michael.hahsler.net/markovDP/reference/act.md),
[`action_state_helpers`](http://michael.hahsler.net/markovDP/reference/action_state_helpers.md),
[`reachable_states()`](http://michael.hahsler.net/markovDP/reference/reachable_states.md),
[`sample_MDP.MDPSample()`](http://michael.hahsler.net/markovDP/reference/sample_MDP.MDPSample.md),
[`start`](http://michael.hahsler.net/markovDP/reference/start.md)

## Examples

``` r
# EXAMPLE
```
