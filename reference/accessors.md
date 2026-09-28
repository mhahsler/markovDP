# Access to Parts of the Model Description

Functions to provide uniform access to different parts of the MDP
problem description.

## Usage

``` r
transition_matrix(
  model,
  action = NULL,
  start.state = NULL,
  end.state = NULL,
  ...,
  sparse = NULL,
  drop = TRUE,
  simplify = FALSE,
  trans_keyword = TRUE
)

reward_matrix(
  model,
  action = NULL,
  start.state = NULL,
  end.state = NULL,
  ...,
  sparse = NULL,
  drop = TRUE,
  simplify = FALSE
)

start_vector(model, start = NULL, sparse = NULL)

normalize_MDP(
  model,
  transition_model = TRUE,
  reward = TRUE,
  start = FALSE,
  sparse = NULL,
  precompute_absorbing = TRUE,
  check_and_fix = FALSE,
  progress = TRUE
)
```

## Arguments

- model:

  A [MDP](http://michael.hahsler.net/markovDP/reference/MDP.md) object.

- action:

  name or index of an action.

- start.state, end.state:

  name or index of the state.

- ...:

  further arguments are passed on.

- sparse:

  logical; use sparse matrix representation? `NULL` decides the
  representation based on the memory it would take to store the faster
  dense representation.

- drop:

  logical; drop matrices to vectors when one row/one column is selected.

- simplify:

  logical; try to simplify action lists into a vector or matrix?

- trans_keyword:

  logical; translate keywords like "uniform" into matrices.

- start:

  logical; convert the start probability distribution into a vector.

- transition_model:

  logical; convert the transition probabilities into a list of matrices.

- reward:

  logical; convert the reward model into a list of matrices.

- precompute_absorbing:

  logical; should absorbing states be precalculated?

- check_and_fix:

  logical; checks the structure of the problem description.

- progress:

  logical; show a progress bar with estimated time for completion.

## Value

A list or a list of lists of matrices.

## Details

Several parts of the MDP description can be defined in different ways.
In particular, the fields `transition_model`, `reward`, and `start` can
be defined using matrices, data frames, keywords, or functions. See
[MDP](http://michael.hahsler.net/markovDP/reference/MDP.md) for details.
The functions provided here give unified access to the data in these
fields and make writing code easier.

### Transition Probabilities \\p(s'\|s,a)\\

`transition_matrix()` accesses the transition model. The complete model
is a list with one element for each action. Each element contains a
states x states matrix with \\s\\ (`start.state`) as rows and \\s'\\
(`end.state`) as columns. Matrices with a low density can be requested
in sparse format (as a
[Matrix::dgRMatrix](https://rdrr.io/pkg/Matrix/man/dgRMatrix-class.html)).
It is recommended to load package `MatrixExtra` to work with sparse
matrices.

### Reward \\r(s,s',a)\\

`reward_matrix()` accesses the reward model. The preferred
representation is a data.frame with the columns `action`, `start.state`,
`end.state`, and `value`. This is a sparse representation.

The dense representation is a list of lists of matrices. The list levels
are \\a\\ (`action`) and \\s\\ (`start.state`). The matrices are column
vectors with rows representing \\s'\\ (`end.state`).

To represent rewards as sparse matrices, **rewards that correspond to a
transition with probability zero are set to zero if the transition
model** is stored as a list of matrices. This makes the reward matrices
as sparse as the transition matrices. The function `normalize_MDP()`
with `sparse = TRUE` will perform this representation.

### Start state

`start_vector()` translates the start state description into a
probability vector.

### Sparse Matrices and Normalizing MDPs

Different components can be specified in various ways. It is often
necessary to convert each component into a specific form (e.g., a dense
matrix) to save time when accessing it. Convert the Complete MDP
Description into a consistent form `normalize_MDP()` converts all
components of the MDP description into a consistent form and returns a
new MDP definition where `transition_model`, `reward`, and `start` are
normalized. This includes the internal representation (dense, sparse, as
a data.frame) and also, `states`, and `actions` are ordered as given in
the problem definition to make safe access using numerical indices
possible. Normalized MDP descriptions can be used in custom code that
expects consistently a certain format.

The default behavior of `sparse = NULL` uses parse matrices for large
models where the dense transition model would need more than
`options("MDP_SPARSE_LIMIT")` (the default is about 100 MB which can be
changed using [`options()`](https://rdrr.io/r/base/options.html)).
Smaller models use faster dense matrices.

## See also

Other MDP:
[`MDP()`](http://michael.hahsler.net/markovDP/reference/MDP.md),
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
[`unreachable_states()`](http://michael.hahsler.net/markovDP/reference/unreachable_states.md)

## Author

Michael Hahsler

## Examples

``` r
data("Maze")
gw_matrix(Maze)
#>      [,1]     [,2]     [,3]     [,4]    
#> [1,] "s(1,1)" "s(1,2)" "s(1,3)" "s(1,4)"
#> [2,] "s(2,1)" NA       "s(2,3)" "s(2,4)"
#> [3,] "s(3,1)" "s(3,2)" "s(3,3)" "s(3,4)"

# here is the internal structure of the Maze object
str(Maze)
#> List of 10
#>  $ name            : chr "Stuart Russell's 3x4 Maze"
#>  $ discount        : num 1
#>  $ horizon         : num Inf
#>  $ states          : chr [1:11] "s(1,1)" "s(2,1)" "s(3,1)" "s(1,2)" ...
#>  $ actions         : chr [1:4] "up" "right" "down" "left"
#>  $ start           : chr "s(3,1)"
#>  $ transition_model:List of 4
#>   ..$ up   : num [1:11, 1:11] 0.9 0.8 0 0.1 0 0 0 0 0 0 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ : chr [1:11] "s(1,1)" "s(2,1)" "s(3,1)" "s(1,2)" ...
#>   .. .. ..$ : chr [1:11] "s(1,1)" "s(2,1)" "s(3,1)" "s(1,2)" ...
#>   ..$ right: num [1:11, 1:11] 0.1 0.1 0 0 0 0 0 0 0 0 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ : chr [1:11] "s(1,1)" "s(2,1)" "s(3,1)" "s(1,2)" ...
#>   .. .. ..$ : chr [1:11] "s(1,1)" "s(2,1)" "s(3,1)" "s(1,2)" ...
#>   ..$ down : num [1:11, 1:11] 0.1 0 0 0.1 0 0 0 0 0 0 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ : chr [1:11] "s(1,1)" "s(2,1)" "s(3,1)" "s(1,2)" ...
#>   .. .. ..$ : chr [1:11] "s(1,1)" "s(2,1)" "s(3,1)" "s(1,2)" ...
#>   ..$ left : num [1:11, 1:11] 0.9 0.1 0 0.8 0 0 0 0 0 0 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ : chr [1:11] "s(1,1)" "s(2,1)" "s(3,1)" "s(1,2)" ...
#>   .. .. ..$ : chr [1:11] "s(1,1)" "s(2,1)" "s(3,1)" "s(1,2)" ...
#>  $ reward          :List of 4
#>   ..$ up   : num [1:11, 1:11] -0.04 -0.04 0 -0.04 0 0 0 0 0 0 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ : chr [1:11] "s(1,1)" "s(2,1)" "s(3,1)" "s(1,2)" ...
#>   .. .. ..$ : chr [1:11] "s(1,1)" "s(2,1)" "s(3,1)" "s(1,2)" ...
#>   ..$ right: num [1:11, 1:11] -0.04 -0.04 0 0 0 0 0 0 0 0 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ : chr [1:11] "s(1,1)" "s(2,1)" "s(3,1)" "s(1,2)" ...
#>   .. .. ..$ : chr [1:11] "s(1,1)" "s(2,1)" "s(3,1)" "s(1,2)" ...
#>   ..$ down : num [1:11, 1:11] -0.04 0 0 -0.04 0 0 0 0 0 0 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ : chr [1:11] "s(1,1)" "s(2,1)" "s(3,1)" "s(1,2)" ...
#>   .. .. ..$ : chr [1:11] "s(1,1)" "s(2,1)" "s(3,1)" "s(1,2)" ...
#>   ..$ left : num [1:11, 1:11] -0.04 -0.04 0 -0.04 0 0 0 0 0 0 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ : chr [1:11] "s(1,1)" "s(2,1)" "s(3,1)" "s(1,2)" ...
#>   .. .. ..$ : chr [1:11] "s(1,1)" "s(2,1)" "s(3,1)" "s(1,2)" ...
#>  $ info            :List of 6
#>   ..$ gridworld       : logi TRUE
#>   ..$ dim             : num [1:2] 3 4
#>   ..$ start           : chr "s(3,1)"
#>   ..$ goal            : chr "s(1,4)"
#>   ..$ state_labels    :List of 3
#>   .. ..$ s(3,1): chr "Start"
#>   .. ..$ s(2,4): chr "-1"
#>   .. ..$ s(1,4): chr "Goal: +1"
#>   ..$ absorbing_states: chr [1:2] "s(1,4)" "s(2,4)"
#>  $ absorbing_states: chr [1:2] "s(1,4)" "s(2,4)"
#>  - attr(*, "class")= chr [1:2] "MDPModel" "MDP"

# List of |A| transition matrices. One per action in the from start.states x end.states
Maze$transition_model
#> $up
#>        s(1,1) s(2,1) s(3,1) s(1,2) s(3,2) s(1,3) s(2,3) s(3,3) s(1,4) s(2,4)
#> s(1,1)    0.9    0.0    0.0    0.1    0.0    0.0    0.0    0.0    0.0    0.0
#> s(2,1)    0.8    0.2    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0
#> s(3,1)    0.0    0.8    0.1    0.0    0.1    0.0    0.0    0.0    0.0    0.0
#> s(1,2)    0.1    0.0    0.0    0.8    0.0    0.1    0.0    0.0    0.0    0.0
#> s(3,2)    0.0    0.0    0.1    0.0    0.8    0.0    0.0    0.1    0.0    0.0
#> s(1,3)    0.0    0.0    0.0    0.1    0.0    0.8    0.0    0.0    0.1    0.0
#> s(2,3)    0.0    0.0    0.0    0.0    0.0    0.8    0.1    0.0    0.0    0.1
#> s(3,3)    0.0    0.0    0.0    0.0    0.1    0.0    0.8    0.0    0.0    0.0
#> s(1,4)    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    1.0    0.0
#> s(2,4)    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    1.0
#> s(3,4)    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.1    0.0    0.8
#>        s(3,4)
#> s(1,1)    0.0
#> s(2,1)    0.0
#> s(3,1)    0.0
#> s(1,2)    0.0
#> s(3,2)    0.0
#> s(1,3)    0.0
#> s(2,3)    0.0
#> s(3,3)    0.1
#> s(1,4)    0.0
#> s(2,4)    0.0
#> s(3,4)    0.1
#> 
#> $right
#>        s(1,1) s(2,1) s(3,1) s(1,2) s(3,2) s(1,3) s(2,3) s(3,3) s(1,4) s(2,4)
#> s(1,1)    0.1    0.1    0.0    0.8    0.0    0.0    0.0    0.0    0.0    0.0
#> s(2,1)    0.1    0.8    0.1    0.0    0.0    0.0    0.0    0.0    0.0    0.0
#> s(3,1)    0.0    0.1    0.1    0.0    0.8    0.0    0.0    0.0    0.0    0.0
#> s(1,2)    0.0    0.0    0.0    0.2    0.0    0.8    0.0    0.0    0.0    0.0
#> s(3,2)    0.0    0.0    0.0    0.0    0.2    0.0    0.0    0.8    0.0    0.0
#> s(1,3)    0.0    0.0    0.0    0.0    0.0    0.1    0.1    0.0    0.8    0.0
#> s(2,3)    0.0    0.0    0.0    0.0    0.0    0.1    0.0    0.1    0.0    0.8
#> s(3,3)    0.0    0.0    0.0    0.0    0.0    0.0    0.1    0.1    0.0    0.0
#> s(1,4)    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    1.0    0.0
#> s(2,4)    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    1.0
#> s(3,4)    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.1
#>        s(3,4)
#> s(1,1)    0.0
#> s(2,1)    0.0
#> s(3,1)    0.0
#> s(1,2)    0.0
#> s(3,2)    0.0
#> s(1,3)    0.0
#> s(2,3)    0.0
#> s(3,3)    0.8
#> s(1,4)    0.0
#> s(2,4)    0.0
#> s(3,4)    0.9
#> 
#> $down
#>        s(1,1) s(2,1) s(3,1) s(1,2) s(3,2) s(1,3) s(2,3) s(3,3) s(1,4) s(2,4)
#> s(1,1)    0.1    0.8    0.0    0.1    0.0    0.0    0.0    0.0    0.0    0.0
#> s(2,1)    0.0    0.2    0.8    0.0    0.0    0.0    0.0    0.0    0.0    0.0
#> s(3,1)    0.0    0.0    0.9    0.0    0.1    0.0    0.0    0.0    0.0    0.0
#> s(1,2)    0.1    0.0    0.0    0.8    0.0    0.1    0.0    0.0    0.0    0.0
#> s(3,2)    0.0    0.0    0.1    0.0    0.8    0.0    0.0    0.1    0.0    0.0
#> s(1,3)    0.0    0.0    0.0    0.1    0.0    0.0    0.8    0.0    0.1    0.0
#> s(2,3)    0.0    0.0    0.0    0.0    0.0    0.0    0.1    0.8    0.0    0.1
#> s(3,3)    0.0    0.0    0.0    0.0    0.1    0.0    0.0    0.8    0.0    0.0
#> s(1,4)    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    1.0    0.0
#> s(2,4)    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    1.0
#> s(3,4)    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.1    0.0    0.0
#>        s(3,4)
#> s(1,1)    0.0
#> s(2,1)    0.0
#> s(3,1)    0.0
#> s(1,2)    0.0
#> s(3,2)    0.0
#> s(1,3)    0.0
#> s(2,3)    0.0
#> s(3,3)    0.1
#> s(1,4)    0.0
#> s(2,4)    0.0
#> s(3,4)    0.9
#> 
#> $left
#>        s(1,1) s(2,1) s(3,1) s(1,2) s(3,2) s(1,3) s(2,3) s(3,3) s(1,4) s(2,4)
#> s(1,1)    0.9    0.1    0.0    0.0    0.0    0.0    0.0    0.0      0    0.0
#> s(2,1)    0.1    0.8    0.1    0.0    0.0    0.0    0.0    0.0      0    0.0
#> s(3,1)    0.0    0.1    0.9    0.0    0.0    0.0    0.0    0.0      0    0.0
#> s(1,2)    0.8    0.0    0.0    0.2    0.0    0.0    0.0    0.0      0    0.0
#> s(3,2)    0.0    0.0    0.8    0.0    0.2    0.0    0.0    0.0      0    0.0
#> s(1,3)    0.0    0.0    0.0    0.8    0.0    0.1    0.1    0.0      0    0.0
#> s(2,3)    0.0    0.0    0.0    0.0    0.0    0.1    0.8    0.1      0    0.0
#> s(3,3)    0.0    0.0    0.0    0.0    0.8    0.0    0.1    0.1      0    0.0
#> s(1,4)    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0      1    0.0
#> s(2,4)    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0      0    1.0
#> s(3,4)    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.8      0    0.1
#>        s(3,4)
#> s(1,1)    0.0
#> s(2,1)    0.0
#> s(3,1)    0.0
#> s(1,2)    0.0
#> s(3,2)    0.0
#> s(1,3)    0.0
#> s(2,3)    0.0
#> s(3,3)    0.0
#> s(1,4)    0.0
#> s(2,4)    0.0
#> s(3,4)    0.1
#> 
transition_matrix(Maze)
#> $up
#>        s(1,1) s(2,1) s(3,1) s(1,2) s(3,2) s(1,3) s(2,3) s(3,3) s(1,4) s(2,4)
#> s(1,1)    0.9    0.0    0.0    0.1    0.0    0.0    0.0    0.0    0.0    0.0
#> s(2,1)    0.8    0.2    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0
#> s(3,1)    0.0    0.8    0.1    0.0    0.1    0.0    0.0    0.0    0.0    0.0
#> s(1,2)    0.1    0.0    0.0    0.8    0.0    0.1    0.0    0.0    0.0    0.0
#> s(3,2)    0.0    0.0    0.1    0.0    0.8    0.0    0.0    0.1    0.0    0.0
#> s(1,3)    0.0    0.0    0.0    0.1    0.0    0.8    0.0    0.0    0.1    0.0
#> s(2,3)    0.0    0.0    0.0    0.0    0.0    0.8    0.1    0.0    0.0    0.1
#> s(3,3)    0.0    0.0    0.0    0.0    0.1    0.0    0.8    0.0    0.0    0.0
#> s(1,4)    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    1.0    0.0
#> s(2,4)    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    1.0
#> s(3,4)    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.1    0.0    0.8
#>        s(3,4)
#> s(1,1)    0.0
#> s(2,1)    0.0
#> s(3,1)    0.0
#> s(1,2)    0.0
#> s(3,2)    0.0
#> s(1,3)    0.0
#> s(2,3)    0.0
#> s(3,3)    0.1
#> s(1,4)    0.0
#> s(2,4)    0.0
#> s(3,4)    0.1
#> 
#> $right
#>        s(1,1) s(2,1) s(3,1) s(1,2) s(3,2) s(1,3) s(2,3) s(3,3) s(1,4) s(2,4)
#> s(1,1)    0.1    0.1    0.0    0.8    0.0    0.0    0.0    0.0    0.0    0.0
#> s(2,1)    0.1    0.8    0.1    0.0    0.0    0.0    0.0    0.0    0.0    0.0
#> s(3,1)    0.0    0.1    0.1    0.0    0.8    0.0    0.0    0.0    0.0    0.0
#> s(1,2)    0.0    0.0    0.0    0.2    0.0    0.8    0.0    0.0    0.0    0.0
#> s(3,2)    0.0    0.0    0.0    0.0    0.2    0.0    0.0    0.8    0.0    0.0
#> s(1,3)    0.0    0.0    0.0    0.0    0.0    0.1    0.1    0.0    0.8    0.0
#> s(2,3)    0.0    0.0    0.0    0.0    0.0    0.1    0.0    0.1    0.0    0.8
#> s(3,3)    0.0    0.0    0.0    0.0    0.0    0.0    0.1    0.1    0.0    0.0
#> s(1,4)    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    1.0    0.0
#> s(2,4)    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    1.0
#> s(3,4)    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.1
#>        s(3,4)
#> s(1,1)    0.0
#> s(2,1)    0.0
#> s(3,1)    0.0
#> s(1,2)    0.0
#> s(3,2)    0.0
#> s(1,3)    0.0
#> s(2,3)    0.0
#> s(3,3)    0.8
#> s(1,4)    0.0
#> s(2,4)    0.0
#> s(3,4)    0.9
#> 
#> $down
#>        s(1,1) s(2,1) s(3,1) s(1,2) s(3,2) s(1,3) s(2,3) s(3,3) s(1,4) s(2,4)
#> s(1,1)    0.1    0.8    0.0    0.1    0.0    0.0    0.0    0.0    0.0    0.0
#> s(2,1)    0.0    0.2    0.8    0.0    0.0    0.0    0.0    0.0    0.0    0.0
#> s(3,1)    0.0    0.0    0.9    0.0    0.1    0.0    0.0    0.0    0.0    0.0
#> s(1,2)    0.1    0.0    0.0    0.8    0.0    0.1    0.0    0.0    0.0    0.0
#> s(3,2)    0.0    0.0    0.1    0.0    0.8    0.0    0.0    0.1    0.0    0.0
#> s(1,3)    0.0    0.0    0.0    0.1    0.0    0.0    0.8    0.0    0.1    0.0
#> s(2,3)    0.0    0.0    0.0    0.0    0.0    0.0    0.1    0.8    0.0    0.1
#> s(3,3)    0.0    0.0    0.0    0.0    0.1    0.0    0.0    0.8    0.0    0.0
#> s(1,4)    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    1.0    0.0
#> s(2,4)    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    1.0
#> s(3,4)    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.1    0.0    0.0
#>        s(3,4)
#> s(1,1)    0.0
#> s(2,1)    0.0
#> s(3,1)    0.0
#> s(1,2)    0.0
#> s(3,2)    0.0
#> s(1,3)    0.0
#> s(2,3)    0.0
#> s(3,3)    0.1
#> s(1,4)    0.0
#> s(2,4)    0.0
#> s(3,4)    0.9
#> 
#> $left
#>        s(1,1) s(2,1) s(3,1) s(1,2) s(3,2) s(1,3) s(2,3) s(3,3) s(1,4) s(2,4)
#> s(1,1)    0.9    0.1    0.0    0.0    0.0    0.0    0.0    0.0      0    0.0
#> s(2,1)    0.1    0.8    0.1    0.0    0.0    0.0    0.0    0.0      0    0.0
#> s(3,1)    0.0    0.1    0.9    0.0    0.0    0.0    0.0    0.0      0    0.0
#> s(1,2)    0.8    0.0    0.0    0.2    0.0    0.0    0.0    0.0      0    0.0
#> s(3,2)    0.0    0.0    0.8    0.0    0.2    0.0    0.0    0.0      0    0.0
#> s(1,3)    0.0    0.0    0.0    0.8    0.0    0.1    0.1    0.0      0    0.0
#> s(2,3)    0.0    0.0    0.0    0.0    0.0    0.1    0.8    0.1      0    0.0
#> s(3,3)    0.0    0.0    0.0    0.0    0.8    0.0    0.1    0.1      0    0.0
#> s(1,4)    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0      1    0.0
#> s(2,4)    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0      0    1.0
#> s(3,4)    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.8      0    0.1
#>        s(3,4)
#> s(1,1)    0.0
#> s(2,1)    0.0
#> s(3,1)    0.0
#> s(1,2)    0.0
#> s(3,2)    0.0
#> s(1,3)    0.0
#> s(2,3)    0.0
#> s(3,3)    0.0
#> s(1,4)    0.0
#> s(2,4)    0.0
#> s(3,4)    0.1
#> 
transition_matrix(Maze, action = "up", sparse = FALSE)
#>        s(1,1) s(2,1) s(3,1) s(1,2) s(3,2) s(1,3) s(2,3) s(3,3) s(1,4) s(2,4)
#> s(1,1)    0.9    0.0    0.0    0.1    0.0    0.0    0.0    0.0    0.0    0.0
#> s(2,1)    0.8    0.2    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0
#> s(3,1)    0.0    0.8    0.1    0.0    0.1    0.0    0.0    0.0    0.0    0.0
#> s(1,2)    0.1    0.0    0.0    0.8    0.0    0.1    0.0    0.0    0.0    0.0
#> s(3,2)    0.0    0.0    0.1    0.0    0.8    0.0    0.0    0.1    0.0    0.0
#> s(1,3)    0.0    0.0    0.0    0.1    0.0    0.8    0.0    0.0    0.1    0.0
#> s(2,3)    0.0    0.0    0.0    0.0    0.0    0.8    0.1    0.0    0.0    0.1
#> s(3,3)    0.0    0.0    0.0    0.0    0.1    0.0    0.8    0.0    0.0    0.0
#> s(1,4)    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    1.0    0.0
#> s(2,4)    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    1.0
#> s(3,4)    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.1    0.0    0.8
#>        s(3,4)
#> s(1,1)    0.0
#> s(2,1)    0.0
#> s(3,1)    0.0
#> s(1,2)    0.0
#> s(3,2)    0.0
#> s(1,3)    0.0
#> s(2,3)    0.0
#> s(3,3)    0.1
#> s(1,4)    0.0
#> s(2,4)    0.0
#> s(3,4)    0.1
transition_matrix(Maze,
  action = "up",
  start.state = "s(3,1)", end.state = "s(2,1)"
)
#> [1] 0.8

# List of list of reward matrices. 1st level is action and second level is the
#  start state in the form of a column vector with elements for end states.
Maze$reward
#> $up
#>        s(1,1) s(2,1) s(3,1) s(1,2) s(3,2) s(1,3) s(2,3) s(3,3) s(1,4) s(2,4)
#> s(1,1)  -0.04   0.00   0.00  -0.04   0.00   0.00   0.00   0.00   0.00   0.00
#> s(2,1)  -0.04  -0.04   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00
#> s(3,1)   0.00  -0.04  -0.04   0.00  -0.04   0.00   0.00   0.00   0.00   0.00
#> s(1,2)  -0.04   0.00   0.00  -0.04   0.00  -0.04   0.00   0.00   0.00   0.00
#> s(3,2)   0.00   0.00  -0.04   0.00  -0.04   0.00   0.00  -0.04   0.00   0.00
#> s(1,3)   0.00   0.00   0.00  -0.04   0.00  -0.04   0.00   0.00   0.96   0.00
#> s(2,3)   0.00   0.00   0.00   0.00   0.00  -0.04  -0.04   0.00   0.00  -1.04
#> s(3,3)   0.00   0.00   0.00   0.00  -0.04   0.00  -0.04   0.00   0.00   0.00
#> s(1,4)   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00
#> s(2,4)   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00
#> s(3,4)   0.00   0.00   0.00   0.00   0.00   0.00   0.00  -0.04   0.00  -1.04
#>        s(3,4)
#> s(1,1)   0.00
#> s(2,1)   0.00
#> s(3,1)   0.00
#> s(1,2)   0.00
#> s(3,2)   0.00
#> s(1,3)   0.00
#> s(2,3)   0.00
#> s(3,3)  -0.04
#> s(1,4)   0.00
#> s(2,4)   0.00
#> s(3,4)  -0.04
#> 
#> $right
#>        s(1,1) s(2,1) s(3,1) s(1,2) s(3,2) s(1,3) s(2,3) s(3,3) s(1,4) s(2,4)
#> s(1,1)  -0.04  -0.04   0.00  -0.04   0.00   0.00   0.00   0.00   0.00   0.00
#> s(2,1)  -0.04  -0.04  -0.04   0.00   0.00   0.00   0.00   0.00   0.00   0.00
#> s(3,1)   0.00  -0.04  -0.04   0.00  -0.04   0.00   0.00   0.00   0.00   0.00
#> s(1,2)   0.00   0.00   0.00  -0.04   0.00  -0.04   0.00   0.00   0.00   0.00
#> s(3,2)   0.00   0.00   0.00   0.00  -0.04   0.00   0.00  -0.04   0.00   0.00
#> s(1,3)   0.00   0.00   0.00   0.00   0.00  -0.04  -0.04   0.00   0.96   0.00
#> s(2,3)   0.00   0.00   0.00   0.00   0.00  -0.04   0.00  -0.04   0.00  -1.04
#> s(3,3)   0.00   0.00   0.00   0.00   0.00   0.00  -0.04  -0.04   0.00   0.00
#> s(1,4)   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00
#> s(2,4)   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00
#> s(3,4)   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00  -1.04
#>        s(3,4)
#> s(1,1)   0.00
#> s(2,1)   0.00
#> s(3,1)   0.00
#> s(1,2)   0.00
#> s(3,2)   0.00
#> s(1,3)   0.00
#> s(2,3)   0.00
#> s(3,3)  -0.04
#> s(1,4)   0.00
#> s(2,4)   0.00
#> s(3,4)  -0.04
#> 
#> $down
#>        s(1,1) s(2,1) s(3,1) s(1,2) s(3,2) s(1,3) s(2,3) s(3,3) s(1,4) s(2,4)
#> s(1,1)  -0.04  -0.04   0.00  -0.04   0.00   0.00   0.00   0.00   0.00   0.00
#> s(2,1)   0.00  -0.04  -0.04   0.00   0.00   0.00   0.00   0.00   0.00   0.00
#> s(3,1)   0.00   0.00  -0.04   0.00  -0.04   0.00   0.00   0.00   0.00   0.00
#> s(1,2)  -0.04   0.00   0.00  -0.04   0.00  -0.04   0.00   0.00   0.00   0.00
#> s(3,2)   0.00   0.00  -0.04   0.00  -0.04   0.00   0.00  -0.04   0.00   0.00
#> s(1,3)   0.00   0.00   0.00  -0.04   0.00   0.00  -0.04   0.00   0.96   0.00
#> s(2,3)   0.00   0.00   0.00   0.00   0.00   0.00  -0.04  -0.04   0.00  -1.04
#> s(3,3)   0.00   0.00   0.00   0.00  -0.04   0.00   0.00  -0.04   0.00   0.00
#> s(1,4)   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00
#> s(2,4)   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00
#> s(3,4)   0.00   0.00   0.00   0.00   0.00   0.00   0.00  -0.04   0.00   0.00
#>        s(3,4)
#> s(1,1)   0.00
#> s(2,1)   0.00
#> s(3,1)   0.00
#> s(1,2)   0.00
#> s(3,2)   0.00
#> s(1,3)   0.00
#> s(2,3)   0.00
#> s(3,3)  -0.04
#> s(1,4)   0.00
#> s(2,4)   0.00
#> s(3,4)  -0.04
#> 
#> $left
#>        s(1,1) s(2,1) s(3,1) s(1,2) s(3,2) s(1,3) s(2,3) s(3,3) s(1,4) s(2,4)
#> s(1,1)  -0.04  -0.04   0.00   0.00   0.00   0.00   0.00   0.00      0   0.00
#> s(2,1)  -0.04  -0.04  -0.04   0.00   0.00   0.00   0.00   0.00      0   0.00
#> s(3,1)   0.00  -0.04  -0.04   0.00   0.00   0.00   0.00   0.00      0   0.00
#> s(1,2)  -0.04   0.00   0.00  -0.04   0.00   0.00   0.00   0.00      0   0.00
#> s(3,2)   0.00   0.00  -0.04   0.00  -0.04   0.00   0.00   0.00      0   0.00
#> s(1,3)   0.00   0.00   0.00  -0.04   0.00  -0.04  -0.04   0.00      0   0.00
#> s(2,3)   0.00   0.00   0.00   0.00   0.00  -0.04  -0.04  -0.04      0   0.00
#> s(3,3)   0.00   0.00   0.00   0.00  -0.04   0.00  -0.04  -0.04      0   0.00
#> s(1,4)   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00      0   0.00
#> s(2,4)   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00      0   0.00
#> s(3,4)   0.00   0.00   0.00   0.00   0.00   0.00   0.00  -0.04      0  -1.04
#>        s(3,4)
#> s(1,1)   0.00
#> s(2,1)   0.00
#> s(3,1)   0.00
#> s(1,2)   0.00
#> s(3,2)   0.00
#> s(1,3)   0.00
#> s(2,3)   0.00
#> s(3,3)   0.00
#> s(1,4)   0.00
#> s(2,4)   0.00
#> s(3,4)  -0.04
#> 
reward_matrix(Maze)
#> $up
#>        s(1,1) s(2,1) s(3,1) s(1,2) s(3,2) s(1,3) s(2,3) s(3,3) s(1,4) s(2,4)
#> s(1,1)  -0.04   0.00   0.00  -0.04   0.00   0.00   0.00   0.00   0.00   0.00
#> s(2,1)  -0.04  -0.04   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00
#> s(3,1)   0.00  -0.04  -0.04   0.00  -0.04   0.00   0.00   0.00   0.00   0.00
#> s(1,2)  -0.04   0.00   0.00  -0.04   0.00  -0.04   0.00   0.00   0.00   0.00
#> s(3,2)   0.00   0.00  -0.04   0.00  -0.04   0.00   0.00  -0.04   0.00   0.00
#> s(1,3)   0.00   0.00   0.00  -0.04   0.00  -0.04   0.00   0.00   0.96   0.00
#> s(2,3)   0.00   0.00   0.00   0.00   0.00  -0.04  -0.04   0.00   0.00  -1.04
#> s(3,3)   0.00   0.00   0.00   0.00  -0.04   0.00  -0.04   0.00   0.00   0.00
#> s(1,4)   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00
#> s(2,4)   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00
#> s(3,4)   0.00   0.00   0.00   0.00   0.00   0.00   0.00  -0.04   0.00  -1.04
#>        s(3,4)
#> s(1,1)   0.00
#> s(2,1)   0.00
#> s(3,1)   0.00
#> s(1,2)   0.00
#> s(3,2)   0.00
#> s(1,3)   0.00
#> s(2,3)   0.00
#> s(3,3)  -0.04
#> s(1,4)   0.00
#> s(2,4)   0.00
#> s(3,4)  -0.04
#> 
#> $right
#>        s(1,1) s(2,1) s(3,1) s(1,2) s(3,2) s(1,3) s(2,3) s(3,3) s(1,4) s(2,4)
#> s(1,1)  -0.04  -0.04   0.00  -0.04   0.00   0.00   0.00   0.00   0.00   0.00
#> s(2,1)  -0.04  -0.04  -0.04   0.00   0.00   0.00   0.00   0.00   0.00   0.00
#> s(3,1)   0.00  -0.04  -0.04   0.00  -0.04   0.00   0.00   0.00   0.00   0.00
#> s(1,2)   0.00   0.00   0.00  -0.04   0.00  -0.04   0.00   0.00   0.00   0.00
#> s(3,2)   0.00   0.00   0.00   0.00  -0.04   0.00   0.00  -0.04   0.00   0.00
#> s(1,3)   0.00   0.00   0.00   0.00   0.00  -0.04  -0.04   0.00   0.96   0.00
#> s(2,3)   0.00   0.00   0.00   0.00   0.00  -0.04   0.00  -0.04   0.00  -1.04
#> s(3,3)   0.00   0.00   0.00   0.00   0.00   0.00  -0.04  -0.04   0.00   0.00
#> s(1,4)   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00
#> s(2,4)   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00
#> s(3,4)   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00  -1.04
#>        s(3,4)
#> s(1,1)   0.00
#> s(2,1)   0.00
#> s(3,1)   0.00
#> s(1,2)   0.00
#> s(3,2)   0.00
#> s(1,3)   0.00
#> s(2,3)   0.00
#> s(3,3)  -0.04
#> s(1,4)   0.00
#> s(2,4)   0.00
#> s(3,4)  -0.04
#> 
#> $down
#>        s(1,1) s(2,1) s(3,1) s(1,2) s(3,2) s(1,3) s(2,3) s(3,3) s(1,4) s(2,4)
#> s(1,1)  -0.04  -0.04   0.00  -0.04   0.00   0.00   0.00   0.00   0.00   0.00
#> s(2,1)   0.00  -0.04  -0.04   0.00   0.00   0.00   0.00   0.00   0.00   0.00
#> s(3,1)   0.00   0.00  -0.04   0.00  -0.04   0.00   0.00   0.00   0.00   0.00
#> s(1,2)  -0.04   0.00   0.00  -0.04   0.00  -0.04   0.00   0.00   0.00   0.00
#> s(3,2)   0.00   0.00  -0.04   0.00  -0.04   0.00   0.00  -0.04   0.00   0.00
#> s(1,3)   0.00   0.00   0.00  -0.04   0.00   0.00  -0.04   0.00   0.96   0.00
#> s(2,3)   0.00   0.00   0.00   0.00   0.00   0.00  -0.04  -0.04   0.00  -1.04
#> s(3,3)   0.00   0.00   0.00   0.00  -0.04   0.00   0.00  -0.04   0.00   0.00
#> s(1,4)   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00
#> s(2,4)   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00
#> s(3,4)   0.00   0.00   0.00   0.00   0.00   0.00   0.00  -0.04   0.00   0.00
#>        s(3,4)
#> s(1,1)   0.00
#> s(2,1)   0.00
#> s(3,1)   0.00
#> s(1,2)   0.00
#> s(3,2)   0.00
#> s(1,3)   0.00
#> s(2,3)   0.00
#> s(3,3)  -0.04
#> s(1,4)   0.00
#> s(2,4)   0.00
#> s(3,4)  -0.04
#> 
#> $left
#>        s(1,1) s(2,1) s(3,1) s(1,2) s(3,2) s(1,3) s(2,3) s(3,3) s(1,4) s(2,4)
#> s(1,1)  -0.04  -0.04   0.00   0.00   0.00   0.00   0.00   0.00      0   0.00
#> s(2,1)  -0.04  -0.04  -0.04   0.00   0.00   0.00   0.00   0.00      0   0.00
#> s(3,1)   0.00  -0.04  -0.04   0.00   0.00   0.00   0.00   0.00      0   0.00
#> s(1,2)  -0.04   0.00   0.00  -0.04   0.00   0.00   0.00   0.00      0   0.00
#> s(3,2)   0.00   0.00  -0.04   0.00  -0.04   0.00   0.00   0.00      0   0.00
#> s(1,3)   0.00   0.00   0.00  -0.04   0.00  -0.04  -0.04   0.00      0   0.00
#> s(2,3)   0.00   0.00   0.00   0.00   0.00  -0.04  -0.04  -0.04      0   0.00
#> s(3,3)   0.00   0.00   0.00   0.00  -0.04   0.00  -0.04  -0.04      0   0.00
#> s(1,4)   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00      0   0.00
#> s(2,4)   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00      0   0.00
#> s(3,4)   0.00   0.00   0.00   0.00   0.00   0.00   0.00  -0.04      0  -1.04
#>        s(3,4)
#> s(1,1)   0.00
#> s(2,1)   0.00
#> s(3,1)   0.00
#> s(1,2)   0.00
#> s(3,2)   0.00
#> s(1,3)   0.00
#> s(2,3)   0.00
#> s(3,3)   0.00
#> s(1,4)   0.00
#> s(2,4)   0.00
#> s(3,4)  -0.04
#> 
reward_matrix(Maze, sparse = TRUE)
#> $up
#> Sparse CSR matrix (class 'dgRMatrix')
#> Dimensions: 11 x 11
#> (25 entries, 20.66% full)
#> 
#> $right
#> Sparse CSR matrix (class 'dgRMatrix')
#> Dimensions: 11 x 11
#> (24 entries, 19.83% full)
#> 
#> $down
#> Sparse CSR matrix (class 'dgRMatrix')
#> Dimensions: 11 x 11
#> (24 entries, 19.83% full)
#> 
#> $left
#> Sparse CSR matrix (class 'dgRMatrix')
#> Dimensions: 11 x 11
#> (23 entries, 19.01% full)
#> 
reward_matrix(Maze,
  action = "up",
  start.state = "s(3,1)", end.state = "s(2,1)"
)
#> [1] -0.04

# Translate the initial start probability vector
Maze$start
#> [1] "s(3,1)"
start_vector(Maze, sparse = FALSE)
#> s(1,1) s(2,1) s(3,1) s(1,2) s(3,2) s(1,3) s(2,3) s(3,3) s(1,4) s(2,4) s(3,4) 
#>      0      0      1      0      0      0      0      0      0      0      0 
start_vector(Maze, sparse = "states")
#> [1] "s(3,1)"
start_vector(Maze, sparse = "index")
#> [1] 3

# Normalize the whole model using sparse representation
Maze_norm <- normalize_MDP(Maze, sparse = TRUE)
str(Maze_norm)
#> List of 10
#>  $ name            : chr "Stuart Russell's 3x4 Maze"
#>  $ discount        : num 1
#>  $ horizon         : num Inf
#>  $ states          : chr [1:11] "s(1,1)" "s(2,1)" "s(3,1)" "s(1,2)" ...
#>  $ actions         : chr [1:4] "up" "right" "down" "left"
#>  $ start           : chr "s(3,1)"
#>  $ transition_model:List of 4
#>   ..$ up   :Formal class 'dgRMatrix' [package "Matrix"] with 6 slots
#>   .. .. ..@ p       : int [1:12] 0 2 4 7 10 13 16 19 22 23 ...
#>   .. .. ..@ j       : int [1:27] 0 3 0 1 1 2 4 0 3 5 ...
#>   .. .. ..@ Dim     : int [1:2] 11 11
#>   .. .. ..@ Dimnames:List of 2
#>   .. .. .. ..$ : NULL
#>   .. .. .. ..$ : NULL
#>   .. .. ..@ x       : num [1:27] 0.9 0.1 0.8 0.2 0.8 0.1 0.1 0.1 0.8 0.1 ...
#>   .. .. ..@ factors : list()
#>   ..$ right:Formal class 'dgRMatrix' [package "Matrix"] with 6 slots
#>   .. .. ..@ p       : int [1:12] 0 3 6 9 11 13 16 19 22 23 ...
#>   .. .. ..@ j       : int [1:26] 0 1 3 0 1 2 1 2 4 3 ...
#>   .. .. ..@ Dim     : int [1:2] 11 11
#>   .. .. ..@ Dimnames:List of 2
#>   .. .. .. ..$ : NULL
#>   .. .. .. ..$ : NULL
#>   .. .. ..@ x       : num [1:26] 0.1 0.1 0.8 0.1 0.8 0.1 0.1 0.1 0.8 0.2 ...
#>   .. .. ..@ factors : list()
#>   ..$ down :Formal class 'dgRMatrix' [package "Matrix"] with 6 slots
#>   .. .. ..@ p       : int [1:12] 0 3 5 7 10 13 16 19 22 23 ...
#>   .. .. ..@ j       : int [1:26] 0 1 3 1 2 2 4 0 3 5 ...
#>   .. .. ..@ Dim     : int [1:2] 11 11
#>   .. .. ..@ Dimnames:List of 2
#>   .. .. .. ..$ : NULL
#>   .. .. .. ..$ : NULL
#>   .. .. ..@ x       : num [1:26] 0.1 0.8 0.1 0.2 0.8 0.9 0.1 0.1 0.8 0.1 ...
#>   .. .. ..@ factors : list()
#>   ..$ left :Formal class 'dgRMatrix' [package "Matrix"] with 6 slots
#>   .. .. ..@ p       : int [1:12] 0 2 5 7 9 11 14 17 20 21 ...
#>   .. .. ..@ j       : int [1:25] 0 1 0 1 2 1 2 0 3 2 ...
#>   .. .. ..@ Dim     : int [1:2] 11 11
#>   .. .. ..@ Dimnames:List of 2
#>   .. .. .. ..$ : NULL
#>   .. .. .. ..$ : NULL
#>   .. .. ..@ x       : num [1:25] 0.9 0.1 0.1 0.8 0.1 0.1 0.9 0.8 0.2 0.8 ...
#>   .. .. ..@ factors : list()
#>  $ reward          :List of 4
#>   ..$ up   :Formal class 'dgRMatrix' [package "Matrix"] with 6 slots
#>   .. .. ..@ p       : int [1:12] 0 2 4 7 10 13 16 19 22 22 ...
#>   .. .. ..@ j       : int [1:25] 0 3 0 1 1 2 4 0 3 5 ...
#>   .. .. ..@ Dim     : int [1:2] 11 11
#>   .. .. ..@ Dimnames:List of 2
#>   .. .. .. ..$ : chr [1:11] "s(1,1)" "s(2,1)" "s(3,1)" "s(1,2)" ...
#>   .. .. .. ..$ : chr [1:11] "s(1,1)" "s(2,1)" "s(3,1)" "s(1,2)" ...
#>   .. .. ..@ x       : num [1:25] -0.04 -0.04 -0.04 -0.04 -0.04 -0.04 -0.04 -0.04 -0.04 -0.04 ...
#>   .. .. ..@ factors : list()
#>   ..$ right:Formal class 'dgRMatrix' [package "Matrix"] with 6 slots
#>   .. .. ..@ p       : int [1:12] 0 3 6 9 11 13 16 19 22 22 ...
#>   .. .. ..@ j       : int [1:24] 0 1 3 0 1 2 1 2 4 3 ...
#>   .. .. ..@ Dim     : int [1:2] 11 11
#>   .. .. ..@ Dimnames:List of 2
#>   .. .. .. ..$ : chr [1:11] "s(1,1)" "s(2,1)" "s(3,1)" "s(1,2)" ...
#>   .. .. .. ..$ : chr [1:11] "s(1,1)" "s(2,1)" "s(3,1)" "s(1,2)" ...
#>   .. .. ..@ x       : num [1:24] -0.04 -0.04 -0.04 -0.04 -0.04 -0.04 -0.04 -0.04 -0.04 -0.04 ...
#>   .. .. ..@ factors : list()
#>   ..$ down :Formal class 'dgRMatrix' [package "Matrix"] with 6 slots
#>   .. .. ..@ p       : int [1:12] 0 3 5 7 10 13 16 19 22 22 ...
#>   .. .. ..@ j       : int [1:24] 0 1 3 1 2 2 4 0 3 5 ...
#>   .. .. ..@ Dim     : int [1:2] 11 11
#>   .. .. ..@ Dimnames:List of 2
#>   .. .. .. ..$ : chr [1:11] "s(1,1)" "s(2,1)" "s(3,1)" "s(1,2)" ...
#>   .. .. .. ..$ : chr [1:11] "s(1,1)" "s(2,1)" "s(3,1)" "s(1,2)" ...
#>   .. .. ..@ x       : num [1:24] -0.04 -0.04 -0.04 -0.04 -0.04 -0.04 -0.04 -0.04 -0.04 -0.04 ...
#>   .. .. ..@ factors : list()
#>   ..$ left :Formal class 'dgRMatrix' [package "Matrix"] with 6 slots
#>   .. .. ..@ p       : int [1:12] 0 2 5 7 9 11 14 17 20 20 ...
#>   .. .. ..@ j       : int [1:23] 0 1 0 1 2 1 2 0 3 2 ...
#>   .. .. ..@ Dim     : int [1:2] 11 11
#>   .. .. ..@ Dimnames:List of 2
#>   .. .. .. ..$ : chr [1:11] "s(1,1)" "s(2,1)" "s(3,1)" "s(1,2)" ...
#>   .. .. .. ..$ : chr [1:11] "s(1,1)" "s(2,1)" "s(3,1)" "s(1,2)" ...
#>   .. .. ..@ x       : num [1:23] -0.04 -0.04 -0.04 -0.04 -0.04 -0.04 -0.04 -0.04 -0.04 -0.04 ...
#>   .. .. ..@ factors : list()
#>  $ info            :List of 6
#>   ..$ gridworld       : logi TRUE
#>   ..$ dim             : num [1:2] 3 4
#>   ..$ start           : chr "s(3,1)"
#>   ..$ goal            : chr "s(1,4)"
#>   ..$ state_labels    :List of 3
#>   .. ..$ s(3,1): chr "Start"
#>   .. ..$ s(2,4): chr "-1"
#>   .. ..$ s(1,4): chr "Goal: +1"
#>   ..$ absorbing_states: chr [1:2] "s(1,4)" "s(2,4)"
#>  $ absorbing_states: chr [1:2] "s(1,4)" "s(2,4)"
#>  - attr(*, "class")= chr [1:2] "MDPModel" "MDP"

# Note to make the reward matrix sparse, all rewards
# for transitions with probability of 0 are zeroed out.
reward_matrix(Maze_norm)
#> $up
#> Sparse CSR matrix (class 'dgRMatrix')
#> Dimensions: 11 x 11
#> (25 entries, 20.66% full)
#> 
#> $right
#> Sparse CSR matrix (class 'dgRMatrix')
#> Dimensions: 11 x 11
#> (24 entries, 19.83% full)
#> 
#> $down
#> Sparse CSR matrix (class 'dgRMatrix')
#> Dimensions: 11 x 11
#> (24 entries, 19.83% full)
#> 
#> $left
#> Sparse CSR matrix (class 'dgRMatrix')
#> Dimensions: 11 x 11
#> (23 entries, 19.01% full)
#> 
```
