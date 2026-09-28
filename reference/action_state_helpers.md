# Conversions for Action and State IDs and Labels

Several helper functions to convert state and action (integer) IDs to
labels and vice versa.

## Usage

``` r
normalize_state(state, model, as = "factor")

normalize_state_id(state, model)

normalize_state_label(state, model)

normalize_state_features(state, model = NULL)

normalize_action(action, model, as = "factor")

normalize_action_id(action, model)

normalize_action_label(action, model)

state2features(state)

features2state(x)

s(...)

get_state_features(model)
```

## Arguments

- state:

  a state in any format

- model:

  an MDP model

- as:

  character; specifies the desired output format

- action:

  an action in any format

- x:

  a state feature vector or a matrix of state feature vectors as rows.

- ...:

  features that should be converted into a row vector used to describe a
  state.

## Value

Functions ending in

- `_factor` return a factor,

- `_id` return an integer id,

- `_label` return a character string,

- `_features` return a state feature matrix,

- no ending return the type specified with parameter `as`.

Other functions:

- `state2features()` returns a feature vector/matrix.

- `features2state(x)` returns a state label in the format
  `s(feature list)`.

- `s()` returns a state features row vector.

## Details

`normalize_state()` and `normalize_action()` convert labels or ids into
a desired standard representation. If only the label or the integer id
(i.e., the index) is needed, the additional functions can be used. These
are typically a lot faster.

To support a factored state representation as feature vectors,
`state2features()`, `feature2states()`, and `get_state_features()` are
provided. `get_state_features()` is only available if the model
explicitly stores a finite state space.

**Note:** A factored state is represented as a **row** vector (matrix
with a single row) for a single state (conveniently created via `s()`)
or a matrix with row vectors for a set of states are used. State labels
are constructed in the form `s(feature1, feature2, ...)`. Factored state
representation is used for value function approximation (see
[`solve_MDP_APPROX()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_APPROX.md))
and for
[MDPSample](http://michael.hahsler.net/markovDP/reference/MDPSample.md)
to describe MDP's via a transition function between factored states.

## See also

Other MDP:
[`MDP()`](http://michael.hahsler.net/markovDP/reference/MDP.md),
[`absorbing_states()`](http://michael.hahsler.net/markovDP/reference/absorbing_states.md),
[`act()`](http://michael.hahsler.net/markovDP/reference/act.md),
[`available_actions()`](http://michael.hahsler.net/markovDP/reference/available_actions.md),
[`find_reachable_states()`](http://michael.hahsler.net/markovDP/reference/find_reachable_states.md),
[`reachable_states()`](http://michael.hahsler.net/markovDP/reference/reachable_states.md),
[`sample_MDP()`](http://michael.hahsler.net/markovDP/reference/sample_MDP.md),
[`sample_MDP.MDPSample()`](http://michael.hahsler.net/markovDP/reference/sample_MDP.MDPSample.md),
[`start`](http://michael.hahsler.net/markovDP/reference/start.md),
[`transition_graph()`](http://michael.hahsler.net/markovDP/reference/transition_graph.md),
[`transition_matrix()`](http://michael.hahsler.net/markovDP/reference/accessors.md),
[`unreachable_states()`](http://michael.hahsler.net/markovDP/reference/unreachable_states.md)

Other MDPSample:
[`MDPSample()`](http://michael.hahsler.net/markovDP/reference/MDPSample.md),
[`absorbing_states()`](http://michael.hahsler.net/markovDP/reference/absorbing_states.md),
[`act()`](http://michael.hahsler.net/markovDP/reference/act.md),
[`reachable_states()`](http://michael.hahsler.net/markovDP/reference/reachable_states.md),
[`sample_MDP.MDPSample()`](http://michael.hahsler.net/markovDP/reference/sample_MDP.MDPSample.md),
[`solve_MDP_PG()`](http://michael.hahsler.net/markovDP/reference/solve_MDP_PG.md),
[`start`](http://michael.hahsler.net/markovDP/reference/start.md)

## Examples

``` r
data(Maze)

# states
normalize_state(1, Maze)
#> [1] s(1,1)
#> 11 Levels: s(1,1) s(2,1) s(3,1) s(1,2) s(3,2) s(1,3) s(2,3) s(3,3) ... s(3,4)
normalize_state(1, Maze, as = "id")
#> [1] 1
normalize_state(1, Maze, as = "label")
#> [1] "s(1,1)"
normalize_state(1, Maze, as = "features")
#>        x1 x2
#> s(1,1)  1  1

get_state_features(Maze)
#>        x1 x2
#> s(1,1)  1  1
#> s(2,1)  2  1
#> s(3,1)  3  1
#> s(1,2)  1  2
#> s(3,2)  3  2
#> s(1,3)  1  3
#> s(2,3)  2  3
#> s(3,3)  3  3
#> s(1,4)  1  4
#> s(2,4)  2  4
#> s(3,4)  3  4

# actions
normalize_action(1, Maze)
#> [1] up
#> Levels: up right down left
normalize_action(1, Maze, as = "id")
#> [1] 1
normalize_action(1, Maze, as = "label")
#> [1] "up"

# state label to feature conversion
state2features("s(1,1)")
#>        x1 x2
#> s(1,1)  1  1
s(1,1)
#>      [,1] [,2]
#> [1,]    1    1
```
