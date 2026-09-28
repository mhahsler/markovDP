# Regret of a Policy and Related Measures

Calculates the regret and related measures for a policy relative to a
benchmark policy.

## Usage

``` r
regret(policy, benchmark, start = NULL, relative = FALSE, ...)

action_discrepancy(
  policy,
  benchmark,
  weighted = FALSE,
  proportion = FALSE,
  states = FALSE
)

value_error(policy, benchmark, type = "RMSVE", weighted = FALSE)
```

## Arguments

- policy:

  a solved MDP containing the policy to calculate the regret for.

- benchmark:

  a solved MDP with the (optimal) policy. Regret is calculated relative
  to this policy.

- start:

  start state distribution. If NULL then the start state of the
  `benchmark` is used.

- relative:

  logical; should the relative regret (regret divided by the reward of
  the benchmark) be calculated?

- ...:

  further arguments are passed on to
  [`expected_return()`](http://michael.hahsler.net/markovDP/reference/expected_return.md).

- weighted:

  logical; should mismatched actions or state value errors be weighted
  by the state visit probability for the benchmark? Rarely or never
  visited states will have now less influence on the measure.

- proportion:

  logical; should the action discrepancy be reported as a proportion of
  states with a different action.

- states:

  logical; return the mismatching state ids.

- type:

  type of error root mean square value error (`"RMSVE"`), mean square
  value error (`"MSVE"`), mean absolute value error (`"MAVE"`), absolute
  value error vector (`"AVE"`), value error vector (`"VE"`).

## Value

- `regret()` returns the regret as a difference of expected long-term
  rewards.

- `action_discrepancy()` returns the number or proportion of diverging
  actions.

- `mean_value_error()` returns the mean squared or absolute difference
  in the value function.

## Details

### Regret

Regret for a policy \\\pi\\ is defined as \$\$v\_\pi(s_0) -
v\_\*(s_0),\$\$ where \\v\_\pi(s_0)\\ represents the expected long-term
state value for following policy \\\pi\\ and the starting in state
\\s_0\\ (or a start distribution). The relative regret is calculated as
\$\$\frac{v\_\pi(s_0) - v\_\*(s_0)}{v\_\*(s_0)}.\$\$

Note that for regret, usually the optimal policy \\\pi^\*\\ is used as
the benchmark. Since the optimal policy may not be known, regret
relative to the best known policy can be used.

Regret is only valid with converged value functions. This means that
either the solver has converged, or the value function was estimated for
the policy using converged
[`policy_evaluation()`](http://michael.hahsler.net/markovDP/reference/policy_evaluation.md).

### Action Discrepancy

The action discrepancy measures the difference between two policies as
the number of states for which the prescribed action in the policies
differs. Often, a policy is compared to the best known policy called the
benchmark policy.

Some times two actions are equivalent (have the same q-value) and the
algorithm breaks the tie randomly. The implementation accounts for this
case.

The action discrepancy can be calculated as a proportion of different
actions or be weighted by the state visit probability given the
benchmark policy. Both weighted and proportional action discrepancy is
scaled in \\\[0, 1\]\\.

### Root Mean Squared Value Error

The root mean value error \$\$\sqrt{\text{VE}} = \sqrt{\|\|v\_\pi -
v\_\*\|\|^2}\$\$ is the sum of the squared differences of state values
between a solution's value function and the optimal value function. For
\\v\_\*\\, the value function of the benchmark solution is used. Related
measures like MSVE (means squared value error), MAVE (means absolute
value error), AVE (absolute value error) and VE (value error) are also
provided.

The error can also be weighted by the state visit probability given the
benchmark policy. This results in the expected error with respect to the
state visit distribution of the benchmark policy. This may be important
to evaluate methods that focus only on estimating the value function for
states that are actually visited using the policy.

## See also

Other policy:
[`action()`](http://michael.hahsler.net/markovDP/reference/action.md),
[`expected_return()`](http://michael.hahsler.net/markovDP/reference/expected_return.md),
[`greedy_action()`](http://michael.hahsler.net/markovDP/reference/greedy_action.md),
[`policy()`](http://michael.hahsler.net/markovDP/reference/policy.md),
[`policy_evaluation()`](http://michael.hahsler.net/markovDP/reference/policy_evaluation.md),
[`visit_probability()`](http://michael.hahsler.net/markovDP/reference/visit_probability.md)

## Author

Michael Hahsler

## Examples

``` r
data(Maze)

sol_optimal <- solve_MDP(Maze)

# a manual policy (go up and in some squares to the right)
acts <- rep("up", times = length(Maze$states))
names(acts) <- Maze$states
acts[c("s(1,1)", "s(1,2)", "s(1,3)")] <- "right"

sol_manual <- add_policy(Maze, manual_policy(Maze, acts, estimate_V = TRUE))

# compare the policies side-by-side
cbind(opt = policy(sol_optimal), manual = policy(sol_manual))
#>    opt.state     opt.V opt.action manual.state   manual.V manual.action
#> 1     s(1,1) 0.8115564      right       s(1,1)  0.8115582         right
#> 2     s(2,1) 0.7615521         up       s(2,1)  0.7615582            up
#> 3     s(3,1) 0.7052527         up       s(3,1)  0.6711367            up
#> 4     s(1,2) 0.8678082      right       s(1,2)  0.8678082         right
#> 5     s(3,2) 0.6551492       left       s(3,2)  0.3488563            up
#> 6     s(1,3) 0.9178082      right       s(1,3)  0.9178082         right
#> 7     s(2,3) 0.6602740         up       s(2,3)  0.6602740            up
#> 8     s(3,3) 0.6110843       left       s(3,3)  0.4345002            up
#> 9     s(1,4) 0.0000000       down       s(1,4)  0.0000000            up
#> 10    s(2,4) 0.0000000      right       s(2,4)  0.0000000            up
#> 11    s(3,4) 0.3872455       left       s(3,4) -0.8850705            up

# the regret is very small. It is about 4.8% of the optimal reward
regret(sol_manual, benchmark = sol_optimal)
#> Warning: model does not contain a converged solution. Using policy evaluation to obtain the value function.
#> [1] 0.03411595
regret(sol_manual, benchmark = sol_optimal, relative = TRUE)
#> Warning: model does not contain a converged solution. Using policy evaluation to obtain the value function.
#> [1] 0.04837408

# The number of different actions (excluding equivalent actions) is 3.
# This about 27% of the actions in the policy. 
action_discrepancy(sol_manual, benchmark = sol_optimal)
#> [1] 3
action_discrepancy(sol_manual, benchmark = sol_optimal, proportion = TRUE)
#> [1] 0.2727273

# Weighted by the probability that a state will be visited shows that
# only 2.3% of the time a different action would be used.
action_discrepancy(sol_manual, benchmark = sol_optimal, weighted = TRUE)
#> [1] 0.0203388

value_error(sol_manual, benchmark = sol_optimal, type = "VE")
#>  [1] 1.843715e-06 6.073922e-06 3.411595e-02 2.768344e-08 3.062929e-01
#>  [6] 6.327144e-09 1.789588e-08 1.765841e-01 0.000000e+00 0.000000e+00
#> [11] 1.272316e+00
value_error(sol_manual, benchmark = sol_optimal, type = "MAVE")
#> [1] 0.1626652
value_error(sol_manual, benchmark = sol_optimal, type = "RMSVE")
#> [1] 0.398286

# Weighting shows that the expected MAVE (expectation taken over the 
# benchmark policy) is rather small.
value_error(sol_manual, benchmark = sol_optimal, type = "MAVE", weighted = TRUE)
#> [1] 0.001071097
```
