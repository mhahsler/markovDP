# Solve MDPs using Monte Carlo Control

Solve MDPs using Monte Carlo control.

## Usage

``` r
solve_MDP_MC(
  model,
  method = "exploring_starts",
  horizon = NULL,
  discount = NULL,
  n = 100,
  Q = NULL,
  epsilon = NULL,
  alpha = NULL,
  first_visit = TRUE,
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

  string; one of the following solution methods:

  - `'exploring_starts'` - on-policy MC control with exploring starts.

  - `'on_policy'` - on-policy MC control with an \\\epsilon\\-greedy
    policy.

  - `'off_policy"'` - off-policy MC control using an \\\epsilon\\-greedy
    behavior policy.

- horizon:

  an integer with the number of epochs for problems with a finite
  planning horizon. If set to `Inf`, the algorithm continues running
  iterations till it converges to the infinite horizon solution. If
  `NULL`, then the horizon specified in `model` will be used.

- discount:

  discount factor in range \\(0, 1\]\\. If `NULL`, then the discount
  factor specified in `model` will be used.

- n:

  number of episodes used for learning.

- Q:

  an initial state-action value matrix. By default an all 0 matrix is
  used.

- epsilon:

  used for the \\\epsilon\\-greedy behavior policies. A scalar value
  between 0 and 1 or a
  [schedule](http://michael.hahsler.net/markovDP/reference/schedule.md).

- alpha:

  step size (learning rate). A scalar value between 0 and 1 or a
  [schedule](http://michael.hahsler.net/markovDP/reference/schedule.md).

- first_visit:

  if `TRUE` then only the first visit of a state/action pair in an
  episode is used to update Q, otherwise, every-visit update is used.

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

## Details

The idea is to estimate the action value function for a policy as the
average of sampled returns.

\$\$q\_\pi(s,a) = \mathbb{E}\_\pi\[R_i\|S_0=s,A_0=a\] \approx
\frac{1}{n} \sum\_{i=1}^n R_i\$\$

Monte Carlo control simulates a whole episode using the current behavior
policy and uses the sampled reward to update the Q values. For on-policy
methods, the behavior policy is updated to be greedy (i.e., optimal)
with respect to the new Q values. Then the next episode is simulated
till the predefined number of episodes is completed.

### Implemented methods

Implemented are the following temporal difference control methods
described in Sutton and Barto (2018).

- **Monte Carlo Control with exploring Starts** learns the optimal
  greedy policy. It uses the same greedy policy for behavior and target
  (on-policy learning). After each episode, the policy is updated to be
  greedy with respect to the current Q values. To make sure all
  states/action pairs are explored, it uses exploring starts meaning
  that new episodes are started at a randomly chosen state using a
  randomly chooses action.

- **On-policy Monte Carlo Control** learns an epsilon-greedy policy
  which it uses for behavior and as the target policy (on-policy
  learning). An epsilon-greedy policy is used to provide exploration.
  For calculating running averages, an update with \\\alpha = 1/n\\ is
  used by default. A different update factor can be set using the
  parameter `alpha` as either a fixed value or a function with the
  signature `function(t, n)` which returns the factor in the range
  \\\[0,1\]\\.

- **Off-policy Monte Carlo Control** uses for behavior an arbitrary soft
  policy (a soft policy has in each state a probability greater than 0
  for all possible actions). We use an epsilon-greedy policy and the
  method learns a greedy policy using importance sampling. Note: This
  method can only learn from the tail of the sampled runs where greedy
  actions are chosen. This means that it is very inefficient in learning
  the beginning portion of long episodes. This problem is especially
  problematic when larger values for \\\epsilon\\ are used.

### Schedules

- epsilon schedule: `t` is increased by each processed episode.

- alpha schedule: `t` is set to the number of times the a Q-value for
  state-action combination was updated.

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
[`solve_MDP_PG()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_PG.md),
[`solve_MDP_SAMP()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_SAMP.md),
[`solve_MDP_TD()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_TD.md)
