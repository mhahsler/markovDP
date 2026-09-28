# Unreachable States

Find or removes unreachable states using the transition model.

## Usage

``` r
unreachable_states(
  model,
  horizon = Inf,
  sparse = "states",
  progress = TRUE,
  ...
)

remove_unreachable_states(model, ...)
```

## Arguments

- model:

  a [MDP](http://michael.hahsler.net/markovDP/reference/MDP.md) object.

- horizon:

  only states that can be reached within the horizon are reachable.

- sparse:

  logical; return a sparse logical vector?

- progress:

  logical; show a progress bar?

- ...:

  further arguments are passed on.

## Value

`unreachable_states()` returns a logical vector indicating the
unreachable states.

`remove_unreachable_states()` returns a model with all unreachable
states removed.

## Details

The function `unreachable_states()` checks if states cannot be reached
from any other state. It performs a depth-first search which can be
slow. The search breaks cycles to avoid an infinite loop. The search
depth can be restricted using `horizon`.

The function `remove_unreachable_states()` simplifies a model by
removing unreachable states from the model description.

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
[`transition_matrix()`](http://michael.hahsler.net/markovDP/reference/accessors.md)

## Author

Michael Hahsler

## Examples

``` r
# create a Maze with an unreachable state

maze_unreach <- gw_read_maze(
    textConnection(c("XXXXXX", 
                     "XS X X",
                     "X  XXX",
                     "X   GX",
                     "XXXXXX")))
gw_plot(maze_unreach)


unreachable_states(maze_unreach)
#> [1] "s(2,5)"
unreachable_states(maze_unreach, sparse = FALSE)
#> s(2,2) s(3,2) s(4,2) s(2,3) s(3,3) s(4,3) s(4,4) s(2,5) s(4,5) 
#>  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE   TRUE  FALSE 

maze <- remove_unreachable_states(maze_unreach)
unreachable_states(maze)
#> character(0)
gw_plot(maze)
```
