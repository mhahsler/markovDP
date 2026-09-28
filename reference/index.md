# Package index

## MDP Models

Define finite-state MDPs and inspect their states, actions, transitions,
and rewards.

- [`MDP()`](http://michael.hahsler.net/markovDP/reference/MDP.md)
  [`S()`](http://michael.hahsler.net/markovDP/reference/MDP.md)
  [`A()`](http://michael.hahsler.net/markovDP/reference/MDP.md)
  [`is_solved_MDP()`](http://michael.hahsler.net/markovDP/reference/MDP.md)
  [`is_converged_MDP()`](http://michael.hahsler.net/markovDP/reference/MDP.md)
  [`P_()`](http://michael.hahsler.net/markovDP/reference/MDP.md)
  [`R_()`](http://michael.hahsler.net/markovDP/reference/MDP.md) :
  Define an MDP Problem with Model Access
- [`MDPSample()`](http://michael.hahsler.net/markovDP/reference/MDPSample.md)
  : Define an MDP With Only Sample Access
- [`available_actions()`](http://michael.hahsler.net/markovDP/reference/available_actions.md)
  : Available Actions in a State
- [`transition_matrix()`](http://michael.hahsler.net/markovDP/reference/accessors.md)
  [`reward_matrix()`](http://michael.hahsler.net/markovDP/reference/accessors.md)
  [`start_vector()`](http://michael.hahsler.net/markovDP/reference/accessors.md)
  [`normalize_MDP()`](http://michael.hahsler.net/markovDP/reference/accessors.md)
  : Access to Parts of the Model Description
- [`act()`](http://michael.hahsler.net/markovDP/reference/act.md) :
  Perform an Action
- [`sample_MDP()`](http://michael.hahsler.net/markovDP/reference/sample_MDP.md)
  : Sample Trajectories from an MDP
- [`absorbing_states()`](http://michael.hahsler.net/markovDP/reference/absorbing_states.md)
  : Absorbing States
- [`normalize_state()`](http://michael.hahsler.net/markovDP/reference/action_state_helpers.md)
  [`normalize_state_id()`](http://michael.hahsler.net/markovDP/reference/action_state_helpers.md)
  [`normalize_state_label()`](http://michael.hahsler.net/markovDP/reference/action_state_helpers.md)
  [`normalize_state_features()`](http://michael.hahsler.net/markovDP/reference/action_state_helpers.md)
  [`normalize_action()`](http://michael.hahsler.net/markovDP/reference/action_state_helpers.md)
  [`normalize_action_id()`](http://michael.hahsler.net/markovDP/reference/action_state_helpers.md)
  [`normalize_action_label()`](http://michael.hahsler.net/markovDP/reference/action_state_helpers.md)
  [`state2features()`](http://michael.hahsler.net/markovDP/reference/action_state_helpers.md)
  [`features2state()`](http://michael.hahsler.net/markovDP/reference/action_state_helpers.md)
  [`s()`](http://michael.hahsler.net/markovDP/reference/action_state_helpers.md)
  [`get_state_features()`](http://michael.hahsler.net/markovDP/reference/action_state_helpers.md)
  : Conversions for Action and State IDs and Labels
- [`find_reachable_states()`](http://michael.hahsler.net/markovDP/reference/find_reachable_states.md)
  : Find Reachable State Space from a Transition Model Function
- [`reachable_states()`](http://michael.hahsler.net/markovDP/reference/reachable_states.md)
  : Find Reachable States
- [`sample_MDP(`*`<MDPSample>`*`)`](http://michael.hahsler.net/markovDP/reference/sample_MDP.MDPSample.md)
  : Sample Trajectories from an MDPSample
- [`start(`*`<MDPModel>`*`)`](http://michael.hahsler.net/markovDP/reference/start.md)
  [`start(`*`<MDPSample>`*`)`](http://michael.hahsler.net/markovDP/reference/start.md)
  : Sample a Start State
- [`transition_graph()`](http://michael.hahsler.net/markovDP/reference/transition_graph.md)
  [`plot_transition_graph()`](http://michael.hahsler.net/markovDP/reference/transition_graph.md)
  [`curve_multiple_directed()`](http://michael.hahsler.net/markovDP/reference/transition_graph.md)
  : Transition Graph
- [`unreachable_states()`](http://michael.hahsler.net/markovDP/reference/unreachable_states.md)
  [`remove_unreachable_states()`](http://michael.hahsler.net/markovDP/reference/unreachable_states.md)
  : Unreachable States

## Solvers

Solve MDPs with dynamic programming, linear programming, sampling, and
reinforcement learning methods.

- [`solve_MDP()`](http://michael.hahsler.net/markovDP/reference/solve_MDP.md)
  : Solve an MDP Problem
- [`solve_MDP_DP()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_DP.md)
  : Solve MDPs using Dynamic Programming
- [`solve_MDP_TD()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_TD.md)
  : Solve MDPs using Tabular Temporal Differencing
- [`solve_MDP_MC()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_MC.md)
  : Solve MDPs using Monte Carlo Control
- [`solve_MDP_LP()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_LP.md)
  : Solve MDPs using Linear Programming
- [`solve_MDP_SAMP()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_SAMP.md)
  : Solve MDPs using Random-Sampling
- [`solve_MDP_APPROX()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_APPROX.md)
  [`approx_Q_value()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_APPROX.md)
  [`approx_greedy_action()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_APPROX.md)
  [`approx_greedy_policy()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_APPROX.md)
  [`approx_V_plot()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_APPROX.md)
  : Solve MDPs with Temporal Differencing with Function Approximation
- [`solve_MDP_PG()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_PG.md)
  : Solve MDPs with Policy Gradient Methods
- [`schedule_exp()`](http://michael.hahsler.net/markovDP/reference/schedule.md)
  [`schedule_exp2()`](http://michael.hahsler.net/markovDP/reference/schedule.md)
  [`schedule_log()`](http://michael.hahsler.net/markovDP/reference/schedule.md)
  [`schedule_linear()`](http://michael.hahsler.net/markovDP/reference/schedule.md)
  [`schedule_harmonic()`](http://michael.hahsler.net/markovDP/reference/schedule.md)
  : Schedules to Reduce Alpha, Epsilon and Other Parameters
- [`convergence_horizon()`](http://michael.hahsler.net/markovDP/reference/convergence_horizon.md)
  : Estimate the Convergence Horizon for an Infinite-Horizon MDP

## Value Functions

Work with value functions.

- [`Q_values()`](http://michael.hahsler.net/markovDP/reference/Q_values.md)
  [`Q_zero()`](http://michael.hahsler.net/markovDP/reference/Q_values.md)
  [`Q_random()`](http://michael.hahsler.net/markovDP/reference/Q_values.md)
  : Q-Values
- [`bellman_update()`](http://michael.hahsler.net/markovDP/reference/bellman_update.md)
  [`bellman_operator()`](http://michael.hahsler.net/markovDP/reference/bellman_update.md)
  : Bellman Update and Bellman operator
- [`value_function()`](http://michael.hahsler.net/markovDP/reference/value_function.md)
  [`plot_value_function()`](http://michael.hahsler.net/markovDP/reference/value_function.md)
  [`V_zero()`](http://michael.hahsler.net/markovDP/reference/value_function.md)
  [`V_random()`](http://michael.hahsler.net/markovDP/reference/value_function.md)
  : Value Function

## Policies

Create and evaluate policies, choose actions, and calculate value
functions and returns.

- [`policy()`](http://michael.hahsler.net/markovDP/reference/policy.md)
  [`add_policy()`](http://michael.hahsler.net/markovDP/reference/policy.md)
  [`random_policy()`](http://michael.hahsler.net/markovDP/reference/policy.md)
  [`manual_policy()`](http://michael.hahsler.net/markovDP/reference/policy.md)
  [`induced_transition_matrix()`](http://michael.hahsler.net/markovDP/reference/policy.md)
  [`induced_reward_matrix()`](http://michael.hahsler.net/markovDP/reference/policy.md)
  : Extract, Create Add a Policy to a Model
- [`action()`](http://michael.hahsler.net/markovDP/reference/action.md)
  : Choose an Action Given a Policy
- [`expected_return()`](http://michael.hahsler.net/markovDP/reference/expected_return.md)
  : Calculate the Expected Return of a Policy
- [`greedy_action()`](http://michael.hahsler.net/markovDP/reference/greedy_action.md)
  [`greedy_policy()`](http://michael.hahsler.net/markovDP/reference/greedy_action.md)
  : Greedy Actions and Policies
- [`policy_evaluation()`](http://michael.hahsler.net/markovDP/reference/policy_evaluation.md)
  [`policy_evaluation_LP()`](http://michael.hahsler.net/markovDP/reference/policy_evaluation.md)
  [`policy_evaluation_MC()`](http://michael.hahsler.net/markovDP/reference/policy_evaluation.md)
  [`policy_evaluation_bellman()`](http://michael.hahsler.net/markovDP/reference/policy_evaluation.md)
  : Policy Evaluation
- [`regret()`](http://michael.hahsler.net/markovDP/reference/regret.md)
  [`action_discrepancy()`](http://michael.hahsler.net/markovDP/reference/regret.md)
  [`value_error()`](http://michael.hahsler.net/markovDP/reference/regret.md)
  : Regret of a Policy and Related Measures
- [`visit_probability()`](http://michael.hahsler.net/markovDP/reference/visit_probability.md)
  : State Visit Probability

## Function Approximation

Approximate value functions and policies with state features and basis
transformations.

- [`solve_MDP_APPROX()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_APPROX.md)
  [`approx_Q_value()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_APPROX.md)
  [`approx_greedy_action()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_APPROX.md)
  [`approx_greedy_policy()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_APPROX.md)
  [`approx_V_plot()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_APPROX.md)
  : Solve MDPs with Temporal Differencing with Function Approximation
- [`q_approx_linear()`](http://michael.hahsler.net/markovDP/reference/linear_function_approximation.md)
  [`v_approx_linear()`](http://michael.hahsler.net/markovDP/reference/linear_function_approximation.md)
  [`pi_approx_linear()`](http://michael.hahsler.net/markovDP/reference/linear_function_approximation.md)
  [`approx_value()`](http://michael.hahsler.net/markovDP/reference/linear_function_approximation.md)
  : Linear Function Approximation
- [`transformation_linear_basis()`](http://michael.hahsler.net/markovDP/reference/transformation.md)
  [`transformation_polynomial_basis()`](http://michael.hahsler.net/markovDP/reference/transformation.md)
  [`transformation_RBF_basis()`](http://michael.hahsler.net/markovDP/reference/transformation.md)
  [`transformation_fourier_basis()`](http://michael.hahsler.net/markovDP/reference/transformation.md)
  [`create_basis_coefs()`](http://michael.hahsler.net/markovDP/reference/transformation.md)
  : Transformation Functions for Linear Function Approximation

## Gridworlds

Create gridworld MDPs and inspect their layouts, transitions, and
solutions.

- [`gw_init()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_s2rc()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_rc2s()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_matrix()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_plot()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_plot_transition_graph()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_animate()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_transition_model()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_transition_model_sparse()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_transition_model_named()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_transition_model_end_state()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_maze_MDP()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_random_maze()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_read_maze()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_path()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  : Helper Functions for Gridworld MDPs
- [`Cliff_walking`](http://michael.hahsler.net/markovDP/reference/Cliff_walking.md)
  [`cliff_walking`](http://michael.hahsler.net/markovDP/reference/Cliff_walking.md)
  : Cliff Walking Gridworld MDP
- [`DynaMaze`](http://michael.hahsler.net/markovDP/reference/DynaMaze.md)
  [`dynamaze`](http://michael.hahsler.net/markovDP/reference/DynaMaze.md)
  : The Dyna Maze
- [`Maze`](http://michael.hahsler.net/markovDP/reference/Maze.md)
  [`maze`](http://michael.hahsler.net/markovDP/reference/Maze.md) :
  Steward Russell's 4x3 Maze Gridworld MDP
- [`Windy_gridworld`](http://michael.hahsler.net/markovDP/reference/Windy_gridworld.md)
  [`windy_gridworld`](http://michael.hahsler.net/markovDP/reference/Windy_gridworld.md)
  : Windy Gridworld MDP Windy Gridworld MDP

## Example MDPs

Explore the included maze, cliff-walking, and windy-gridworld examples.

- [`Maze`](http://michael.hahsler.net/markovDP/reference/Maze.md)
  [`maze`](http://michael.hahsler.net/markovDP/reference/Maze.md) :
  Steward Russell's 4x3 Maze Gridworld MDP
- [`Cliff_walking`](http://michael.hahsler.net/markovDP/reference/Cliff_walking.md)
  [`cliff_walking`](http://michael.hahsler.net/markovDP/reference/Cliff_walking.md)
  : Cliff Walking Gridworld MDP
- [`Windy_gridworld`](http://michael.hahsler.net/markovDP/reference/Windy_gridworld.md)
  [`windy_gridworld`](http://michael.hahsler.net/markovDP/reference/Windy_gridworld.md)
  : Windy Gridworld MDP Windy Gridworld MDP
- [`DynaMaze`](http://michael.hahsler.net/markovDP/reference/DynaMaze.md)
  [`dynamaze`](http://michael.hahsler.net/markovDP/reference/DynaMaze.md)
  : The Dyna Maze

## Visualization

Plot gridworlds and transition graphs using the package’s color
palettes.

- [`gw_init()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_s2rc()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_rc2s()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_matrix()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_plot()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_plot_transition_graph()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_animate()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_transition_model()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_transition_model_sparse()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_transition_model_named()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_transition_model_end_state()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_maze_MDP()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_random_maze()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_read_maze()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  [`gw_path()`](http://michael.hahsler.net/markovDP/reference/gridworld.md)
  : Helper Functions for Gridworld MDPs
- [`transition_graph()`](http://michael.hahsler.net/markovDP/reference/transition_graph.md)
  [`plot_transition_graph()`](http://michael.hahsler.net/markovDP/reference/transition_graph.md)
  [`curve_multiple_directed()`](http://michael.hahsler.net/markovDP/reference/transition_graph.md)
  : Transition Graph

## Utilities

Helper functions.

- [`colors_discrete()`](http://michael.hahsler.net/markovDP/reference/colors.md)
  [`colors_continuous()`](http://michael.hahsler.net/markovDP/reference/colors.md)
  : Default Colors for Visualization
- [`round_stochastic()`](http://michael.hahsler.net/markovDP/reference/round_stochastic.md)
  : Round a stochastic vector or a row-stochastic matrix
