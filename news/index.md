# Changelog

## markovDP 0.99.0 (unpublished)

- Updated the test suite to testthat edition 3 and added tests for the
  main model, solver, policy, sampling, and gridworld workflows.
- Standardized policy data frames on the documented `state` column and
  removed partial-match warnings in solvers, plotting, and sampling.
- Corrected the `transition_model` argument name in
  [`find_reachable_states()`](http://michael.hahsler.net/markovDP/reference/find_reachable_states.md).
- Added RL algorithms
- Renamed the policy return calculation to
  [`expected_return()`](http://michael.hahsler.net/markovDP/reference/expected_return.md)
  and the simulation summary field to `avg_return`.
- Separated code from package pomdp (10/01/2024).
