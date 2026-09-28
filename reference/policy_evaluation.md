# Policy Evaluation

Estimate the value function for a policy.

## Usage

``` r
policy_evaluation(
  model,
  pi = NULL,
  method = "bellman",
  ...,
  progress = TRUE,
  verbose = FALSE
)

policy_evaluation_LP(
  model,
  pi = NULL,
  inf = 1000,
  ...,
  progress = FALSE,
  verbose = FALSE
)

policy_evaluation_MC(
  model,
  pi = NULL,
  n = 1000,
  horizon = NULL,
  first_visit = TRUE,
  ...,
  progress = TRUE,
  verbose = FALSE
)

policy_evaluation_bellman(
  model,
  pi = NULL,
  V = NULL,
  k_backups = 1000L,
  theta = 0.001,
  progress = TRUE,
  verbose = FALSE
)
```

## Arguments

- model:

  an MDP problem specification.

- pi:

  a policy as a data.frame with at least columns for states and action.
  If `NULL`, then the policy in model is used.

- method:

  used method: `"bellman"`, `"LP"`, or `"MC"`.

- ...:

  further parameters are passed on to
  [`sample_MDP()`](http://michael.hahsler.net/markovDP/reference/sample_MDP.md).

- progress:

  logical; show a progress bar with estimated time for completion.

- verbose:

  logical; should progress and approximation errors be printed.

- inf:

  value used to replace infinity for
  [`lpSolve::lp()`](https://rdrr.io/pkg/lpSolve/man/lp.html).

- n:

  number of simulated episodes.

- horizon:

  maximum horizon for episodes.

- first_visit:

  if `TRUE` then only the first visit of a state/action pair in an
  episode is used to update Q, otherwise, every-visit update is used.

- V:

  a vector with estimated state values representing a value function. If
  `model` is a solved model, then the state values are taken from the
  solution.

- k_backups:

  number of look ahead steps used for approximate policy evaluation used
  by the policy iteration method. Set k_backups to `Inf` to only use
  \\\theta\\ as the stopping criterion.

- theta:

  stop when the largest state Bellman error (\\\delta = V\_{k+1} - V\\)
  is less than \\\theta\\.

## Value

a vector with (approximate) state values (U).

## Details

The LP implementation only works for infinite-horizon problems with a
discount factor \<1. It solves the following LP:

\$\$\min \sum\_{s\in S} v(s)\$\$ s.t. \$\$v(s) \ge \sum\_{s' \in S} p(s'
\| s, \pi(s)) \[r(s, \pi(s), s') + \gamma v(s')\],\\ \forall s \in S
\$\$

## Policy Evaluation using Bellman Operator

The value function for a policy can be estimated (called policy
evaluation) by repeatedly applying the Bellman operator \$\$v \leftarrow
B\_\pi(v)\$\$ till convergence.

In each iteration, all state values are updated. In this implementation
updating is stopped when the largest state Bellman error is below a
threshold.

\$\$\|\|v\_{k+1} - v_k\|\|\_\infty \< \theta.\$\$

Or if `k_backups` iterations have been completed.

## References

Sutton, R. S., Barto, A. G. (2020). Reinforcement Learning: An
Introduction. Second edition. The MIT Press.

## See also

Other policy:
[`action()`](http://michael.hahsler.net/markovDP/reference/action.md),
[`expected_return()`](http://michael.hahsler.net/markovDP/reference/expected_return.md),
[`greedy_action()`](http://michael.hahsler.net/markovDP/reference/greedy_action.md),
[`policy()`](http://michael.hahsler.net/markovDP/reference/policy.md),
[`regret()`](http://michael.hahsler.net/markovDP/reference/regret.md),
[`visit_probability()`](http://michael.hahsler.net/markovDP/reference/visit_probability.md)

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

# create several policies:
# 1. optimal policy using value iteration
maze_solved <- solve_MDP(Maze, method = "DP:VI")
pi_opt <- policy(maze_solved)
pi_opt
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

# 2. a manual policy (go up and in some squares to the right)
acts <- rep("up", times = length(Maze$states))
names(acts) <- Maze$states
acts[c("s(1,1)", "s(1,2)", "s(1,3)")] <- "right"
pi_manual <- manual_policy(Maze, acts)
pi_manual
#>     state  V action
#> 1  s(1,1) NA  right
#> 2  s(2,1) NA     up
#> 3  s(3,1) NA     up
#> 4  s(1,2) NA  right
#> 5  s(3,2) NA     up
#> 6  s(1,3) NA  right
#> 7  s(2,3) NA     up
#> 8  s(3,3) NA     up
#> 9  s(1,4) NA     up
#> 10 s(2,4) NA     up
#> 11 s(3,4) NA     up

# 3. a random policy
set.seed(1234)
pi_random <- random_policy(Maze, prob = c(up = .7, right = .1, down = .1, left = 0.1))
pi_random
#>     state  V action
#> 1  s(1,1) NA     up
#> 2  s(2,1) NA     up
#> 3  s(3,1) NA     up
#> 4  s(1,2) NA     up
#> 5  s(3,2) NA   left
#> 6  s(1,3) NA     up
#> 7  s(2,3) NA     up
#> 8  s(3,3) NA     up
#> 9  s(1,4) NA     up
#> 10 s(2,4) NA     up
#> 11 s(3,4) NA     up

# 4. an improved policy based on one policy evaluation and
#   policy improvement step.
V <- policy_evaluation(Maze, pi_random)
Q <- Q_values(Maze, V)
pi_greedy <- greedy_policy(Q)
pi_greedy
#>     state          V action
#> 1  s(1,1) -1.0870605  right
#> 2  s(2,1) -1.4048245     up
#> 3  s(3,1) -1.4597764     up
#> 4  s(1,2) -0.3767883  right
#> 5  s(3,2) -0.7749966  right
#> 6  s(1,3)  0.7104851  right
#> 7  s(2,3) -0.3155125     up
#> 8  s(3,3) -0.5422340     up
#> 9  s(1,4)  0.0000000   down
#> 10 s(2,4)  0.0000000  right
#> 11 s(3,4) -0.6728186   left

# compare the approx. value functions for the policies (we restrict
#    the number of backups for the random policy since it may not converge)
rbind(
  random = policy_evaluation(Maze, pi_random, k_backups = 100),
  manual = policy_evaluation(Maze, pi_manual),
  greedy = policy_evaluation(Maze, pi_greedy),
  optimal = policy_evaluation(Maze, pi_opt)
)
#>             s(1,1)     s(2,1)     s(3,1)     s(1,2)     s(3,2)     s(1,3)
#> random  -1.2316877 -1.2774149 -1.3287059 -0.8650240 -1.3741987 -0.1250940
#> manual   0.8115582  0.7615582  0.6711367  0.8678082  0.3488563  0.9178082
#> greedy   0.8115080  0.7613962  0.6905663  0.8678068  0.5263053  0.9178079
#> optimal  0.8115563  0.7615518  0.7052469  0.8678082  0.6551312  0.9178082
#>             s(2,3)     s(3,3) s(1,4) s(2,4)     s(3,4)
#> random  -0.2652520 -0.4868600      0      0 -0.9872413
#> manual   0.6602740  0.4345002      0      0 -0.8850705
#> greedy   0.6602731  0.5764753      0      0  0.3567542
#> optimal  0.6602740  0.6110441      0      0  0.3871554

# use fist-visit Monte Carlo prediction with 100 episodes 
#   and a max horizon of 100
policy_evaluation(Maze, pi_opt, method = "MC", n = 100, horizon = 100)
#>    s(1,1)    s(2,1)    s(3,1)    s(1,2)    s(3,2)    s(1,3)    s(2,3)    s(3,3) 
#> 0.8387402 0.7845714 0.7336508 0.8950376 0.7088889 0.9459016 0.9040000 0.0000000 
#>    s(1,4)    s(2,4)    s(3,4) 
#> 0.0000000 0.0000000 0.0000000 
```
