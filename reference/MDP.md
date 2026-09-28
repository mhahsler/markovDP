# Define an MDP Problem with Model Access

Defines all the elements of a discrete-time finite state-space MDP
problem.

## Usage

``` r
MDP(
  states,
  actions,
  transition_model,
  reward,
  discount = 0.9,
  horizon = Inf,
  start = "uniform",
  info = NULL,
  name = NA
)

S(model)

A(model)

is_solved_MDP(model, policy = TRUE, approx = FALSE, stop = FALSE)

is_converged_MDP(model, stop = FALSE)

P_(action = NA, start.state = NA, end.state = NA, probability)

R_(action = NA, start.state = NA, end.state = NA, value)
```

## Arguments

- states:

  a character vector specifying the names of the states.

- actions:

  a character vector specifying the names of the available actions.

- transition_model:

  Specifies the transition probabilities between states.

- reward:

  Specifies the rewards, which may depend on the action and state.

- discount:

  numeric; discount rate between 0 and 1.

- horizon:

  numeric; Number of epochs. `Inf` specifies an infinite horizon.

- start:

  Specifies in which state the MDP starts.

- info:

  A list with additional information.

- name:

  a string to identify the MDP problem.

- model:

  an `MDP` object.

- policy:

  logical; solution is an explicit policy.

- approx:

  logical; solution is an approximation function.

- stop:

  logical; stop with an error.

- action:

  an action label or integer. The value `NA` matches any action.

- start.state, end.state:

  state as a state label or an integer. The value `NA` matches any
  state.

- probability, value:

  Values used in the helper functions `P_()` and `R_()`.

## Value

The function returns an object of class MDP which is list with the model
specification.
[`solve_MDP()`](http://michael.hahsler.net/markovDP/reference/solve_MDP.md)
reads the object and adds a list element called `'solution'`.

## Details

Markov decision processes (MDPs) are discrete-time stochastic control
processes. This package implements MDPs with a finite state space.
`MDP()` defines all the elements of an MDP problem, including the
discount rate, the set of states, the set of actions, the transition
model, and the reward model.

We use the following notation. An MDP is a five-tuple:

\\(S,A,P,R, \gamma)\\.

\\S\\ is the set of states; \\A\\ is the set of actions; \\P\\ is the
conditional transition probability matrix between states; \\R\\ is the
reward function; and \\\gamma\\ is the discount factor. We will use
lower case letters to represent a member of a set, e.g., \\s\\ is a
specific state. To refer to the size of a set we will use cardinality,
e.g., the number of actions is \\\|A\|\\.

### Names used for mathematical symbols in code

- \\S, s, s'\\: `'states', start.state', 'end.state'`

- \\A, a\\: `'actions', 'action'`

State names and actions can be specified as strings or index numbers
(e.g., `start.state` can be specified as the index of the state in
`states`). For the specification as data.frames below, `NA` can be used
to mean any `start.state`, `end.state` or `action`.

### Specification of transition model: \\P(s' \| s, a)\\

Transition probability to transition to state \\s'\\ from given state
\\s\\ and action \\a\\. The transition probabilities can be specified in
the following ways:

- A data.frame with columns exactly like the arguments of `P_()`. You
  can use [`rbind()`](https://rdrr.io/r/base/cbind.html) with helper
  function `P_()` to create this data frame. Probabilities can be
  specified multiple times and the definition that appears last in the
  data.frame will take affect.

- A named list of matrices, one for each action. Each matrix is square
  with rows representing start states \\s\\ and columns representing end
  states \\s'\\. Instead of a matrix, also the strings `'identity'` or
  `'uniform'` can be specified.

- A function with the following arguments:

  - A function with the argument list `model`, `action`, `start.state`,
    `end.state` which returns a single transition probability.

  - A function with the argument list `model`, `action`, `start.state`
    which returns a transition probability vector for all end states.
    This vector can be dense, a
    [Matrix::sparseVector](https://rdrr.io/pkg/Matrix/man/sparseVector.html)
    or a named vector only containing the non-zero probabilities named
    by the corresponding end state.

  The arguments `action`, `start.state`, and `end.state` will be always
  called with the state names as a character vectors of length 1.

### Specification of the reward function: \\R(a, s, s')\\

The reward function can be specified in the following ways:

- A data frame with columns named exactly like the arguments of `R_()`.
  You can use [`rbind()`](https://rdrr.io/r/base/cbind.html) with helper
  function `R_()` to create this data frame. Rewards can be specified
  multiple times and the definition that appears last in the data.frame
  will take affect.

- A named list of matrices, one for each action. Each matrix is square
  with rows representing start states \\s\\ and columns representing end
  states \\s'\\.

- A function following the same rules as for transition probabilities.

To avoid overflow problems with rewards, reward values should stay well
within the range of `[-1e10, +1e10]`. `-Inf` can be used as the reward
for unavailable actions and will be translated into a large negative
reward for solvers that only support finite reward values.

### Specification of the Start State

The start state of the agent can be a single state or a distribution
over the states. The start state definition is used as the default when
the return is calculated by
[`expected_return()`](http://michael.hahsler.net/markovDP/reference/expected_return.md)
and for sampling with
[`sample_MDP()`](http://michael.hahsler.net/markovDP/reference/sample_MDP.md).

Options to specify the start state are:

- A string specifying the name of a single starting state.

- An integer in the range \\1\\ to \\n\\ to specify the index of a
  single starting state.

- The string `"uniform"` where the start state is chosen using a uniform
  distribution over all states.

- A probability distribution over the states. That is, a vector of
  \\\|S\|\\ probabilities, that add up to \\1\\.

By default, the start state is selected uniformly from all states.

### Accessing Elements of the MDP

The convenience functions `S()` and `A()` return the set of states and
actions.

See
[accessors](http://michael.hahsler.net/markovDP/reference/accessors.md)
for accessing transition probabilities, rewards, and the start state
distribution.

## See also

Other MDP:
[`absorbing_states()`](http://michael.hahsler.net/markovDP/reference/absorbing_states.md),
[`act()`](http://michael.hahsler.net/markovDP/reference/act.md),
[`action_state_helpers`](http://michael.hahsler.net/markovDP/reference/action_state_helpers.md),
[`available_actions()`](http://michael.hahsler.net/markovDP/reference/available_actions.md),
[`find_reachable_states()`](http://michael.hahsler.net/markovDP/reference/find_reachable_states.md),
[`reachable_states()`](http://michael.hahsler.net/markovDP/reference/reachable_states.md),
[`sample_MDP()`](http://michael.hahsler.net/markovDP/reference/sample_MDP.md),
[`sample_MDP.MDPSample()`](http://michael.hahsler.net/markovDP/reference/sample_MDP.MDPSample.md),
[`start`](http://michael.hahsler.net/markovDP/reference/start.md),
[`transition_graph()`](http://michael.hahsler.net/markovDP/reference/transition_graph.md),
[`transition_matrix()`](http://michael.hahsler.net/markovDP/reference/accessors.md),
[`unreachable_states()`](http://michael.hahsler.net/markovDP/reference/unreachable_states.md)

Other MDP_examples:
[`Cliff_walking`](http://michael.hahsler.net/markovDP/reference/Cliff_walking.md),
[`DynaMaze`](http://michael.hahsler.net/markovDP/reference/DynaMaze.md),
[`Maze`](http://michael.hahsler.net/markovDP/reference/Maze.md),
[`Windy_gridworld`](http://michael.hahsler.net/markovDP/reference/Windy_gridworld.md)

## Author

Michael Hahsler

## Examples

``` r
# simple MDP example
#
# states:    s1 s2 s3 s4
# transitions: forward moves -> and backward moves <-
# start: s1
# reward: s1, s2, s4 = 0 and s3 = 1

car <- MDP(
  states = c("s1", "s2", "s3", "s4"),
  actions = c("forward", "back", "stop"),
  transition <- list(
    forward = rbind(c(0, 1, 0, 0), 
                    c(0, 0, 1, 0), 
                    c(0, 0, 0, 1), 
                    c(0, 0, 0, 1)),
    back =    rbind(c(1, 0, 0, 0), 
                    c(1, 0, 0, 0), 
                    c(0, 1, 0, 0), 
                    c(0, 0, 1, 0)),
    stop = "identity"
  ),
  reward = rbind(
    R_(value = 0),
    R_(end.state = "s3", value = 1)
  ),
  discount = 0.9,
  start = "s1",
  name = "Simple Car MDP"
)

car
#> MDPModel, MDP - Simple Car MDP
#>   Discount factor: 0.9
#>   Horizon: Inf epochs
#>   Size: 3 actions / 4 states
#>   Storage: transition prob as matrix / reward as data.frame. Total size: 6.1 Kb
#>   Start: s1
#>   Model list components: ‘name’, ‘discount’, ‘horizon’, ‘states’,
#>     ‘actions’, ‘start’, ‘transition_model’, ‘reward’, ‘info’

# internal representation
str(car)
#> List of 9
#>  $ name            : chr "Simple Car MDP"
#>  $ discount        : num 0.9
#>  $ horizon         : num Inf
#>  $ states          : chr [1:4] "s1" "s2" "s3" "s4"
#>  $ actions         : chr [1:3] "forward" "back" "stop"
#>  $ start           : chr "s1"
#>  $ transition_model:List of 3
#>   ..$ forward: num [1:4, 1:4] 0 0 0 0 1 0 0 0 0 1 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ : chr [1:4] "s1" "s2" "s3" "s4"
#>   .. .. ..$ : chr [1:4] "s1" "s2" "s3" "s4"
#>   ..$ back   : num [1:4, 1:4] 1 1 0 0 0 0 1 0 0 0 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ : chr [1:4] "s1" "s2" "s3" "s4"
#>   .. .. ..$ : chr [1:4] "s1" "s2" "s3" "s4"
#>   ..$ stop   : chr "identity"
#>  $ reward          :'data.frame':    2 obs. of  4 variables:
#>   ..$ action     : chr [1:2] NA NA
#>   ..$ start.state: chr [1:2] NA NA
#>   ..$ end.state  : chr [1:2] NA "s3"
#>   ..$ value      : int [1:2] 0 1
#>  $ info            : NULL
#>  - attr(*, "class")= chr [1:2] "MDPModel" "MDP"

# accessing elements
S(car)
#> [1] "s1" "s2" "s3" "s4"
A(car)
#> [1] "forward" "back"    "stop"   
start_vector(car, sparse = "states")
#> [1] "s1"
transition_matrix(car)
#> $forward
#>    s1 s2 s3 s4
#> s1  0  1  0  0
#> s2  0  0  1  0
#> s3  0  0  0  1
#> s4  0  0  0  1
#> 
#> $back
#>    s1 s2 s3 s4
#> s1  1  0  0  0
#> s2  1  0  0  0
#> s3  0  1  0  0
#> s4  0  0  1  0
#> 
#> $stop
#> Sparse CSR matrix (class 'dgRMatrix')
#> Dimensions: 4 x 4
#> (4 entries, 25.00% full)
#> 
transition_matrix(car, sparse = TRUE)
#> $forward
#> Sparse CSR matrix (class 'dgRMatrix')
#> Dimensions: 4 x 4
#> (4 entries, 25.00% full)
#> 
#> $back
#> Sparse CSR matrix (class 'dgRMatrix')
#> Dimensions: 4 x 4
#> (4 entries, 25.00% full)
#> 
#> $stop
#> Sparse CSR matrix (class 'dgRMatrix')
#> Dimensions: 4 x 4
#> (4 entries, 25.00% full)
#> 
reward_matrix(car)
#> $forward
#>    s1 s2 s3 s4
#> s1  0  0  0  0
#> s2  0  0  1  0
#> s3  0  0  0  0
#> s4  0  0  0  0
#> 
#> $back
#>    s1 s2 s3 s4
#> s1  0  0  0  0
#> s2  0  0  0  0
#> s3  0  0  0  0
#> s4  0  0  1  0
#> 
#> $stop
#>    s1 s2 s3 s4
#> s1  0  0  1  0
#> s2  0  0  1  0
#> s3  0  0  1  0
#> s4  0  0  1  0
#> 
reward_matrix(car, sparse = TRUE)
#> $forward
#> Sparse CSR matrix (class 'dgRMatrix')
#> Dimensions: 4 x 4
#> (1 entries, 6.25% full)
#> 
#> $back
#> Sparse CSR matrix (class 'dgRMatrix')
#> Dimensions: 4 x 4
#> (1 entries, 6.25% full)
#> 
#> $stop
#> Sparse CSR matrix (class 'dgRMatrix')
#> Dimensions: 4 x 4
#> (4 entries, 25.00% full)
#> 

sol <- solve_MDP(car)
sol
#> MDPModel, MDP - Simple Car MDP
#>   Discount factor: 0.9
#>   Horizon: Inf epochs
#>   Size: 3 actions / 4 states
#>   Storage: transition prob as matrix / reward as matrix. Total size: 12.4 Kb
#>   Start: s1
#>   Model list components: ‘name’, ‘discount’, ‘horizon’, ‘states’,
#>     ‘actions’, ‘start’, ‘transition_model’, ‘reward’, ‘info’,
#>     ‘solution’
#> 
#>   Solved:
#>     Method: ‘VI’
#>     Solution converged: TRUE
#>   Solution list components: ‘method’, ‘policy’, ‘converged’, ‘delta’,
#>     ‘iterations’

policy(sol)
#>   state       V  action
#> 1    s1 8.99906 forward
#> 2    s2 9.99906 forward
#> 3    s3 9.99906    stop
#> 4    s4 9.99906    back
```
