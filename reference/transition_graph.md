# Transition Graph

Returns the transition model as an igraph object.

## Usage

``` r
transition_graph(
  x,
  action = NULL,
  state_col = NULL,
  simplify_transitions = TRUE,
  remove_unavailable_actions = TRUE
)

plot_transition_graph(
  x,
  action = NULL,
  state_col = NULL,
  simplify_transitions = TRUE,
  main = NULL,
  ...
)

curve_multiple_directed(graph, start = 0.3)
```

## Arguments

- x:

  object of class
  [MDP](http://michael.hahsler.net/markovDP/reference/MDP.md).

- action:

  the name or id of an action or a set of actions. By default the
  transition model for all actions is returned.

- state_col:

  colors used to represent the states.

- simplify_transitions:

  logical; combine parallel transition arcs into a single arc.

- remove_unavailable_actions:

  logical; don't show arrows for unavailable actions.

- main:

  a main title for the plot.

- ...:

  further arguments are passed on to
  [`igraph::plot.igraph()`](https://r.igraph.org/reference/plot.igraph.html).

- graph:

  The input graph.

- start:

  The curvature at the two extreme edges.

## Value

returns the transition model as an igraph object.

## Details

The transition model of an MDP is a Markov chain. This function extracts
it as an igraph object.

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
[`transition_matrix()`](http://michael.hahsler.net/markovDP/reference/accessors.md),
[`unreachable_states()`](http://michael.hahsler.net/markovDP/reference/unreachable_states.md)

Other visualization:
[`gridworld`](http://michael.hahsler.net/markovDP/reference/gridworld.md)

## Examples

``` r
data("Maze")

g <- transition_graph(Maze)
g
#> IGRAPH c8e07ef DN-- 11 32 -- 
#> + attr: name (v/c), color (v/c), label (e/c)
#> + edges from c8e07ef (vertex names):
#>  [1] s(1,1)->s(1,1) s(1,1)->s(2,1) s(1,1)->s(1,2) s(2,1)->s(1,1) s(2,1)->s(2,1)
#>  [6] s(2,1)->s(3,1) s(3,1)->s(2,1) s(3,1)->s(3,1) s(3,1)->s(3,2) s(1,2)->s(1,1)
#> [11] s(1,2)->s(1,2) s(1,2)->s(1,3) s(3,2)->s(3,1) s(3,2)->s(3,2) s(3,2)->s(3,3)
#> [16] s(1,3)->s(1,2) s(1,3)->s(1,3) s(1,3)->s(2,3) s(1,3)->s(1,4) s(2,3)->s(1,3)
#> [21] s(2,3)->s(2,3) s(2,3)->s(3,3) s(2,3)->s(2,4) s(3,3)->s(3,2) s(3,3)->s(2,3)
#> [26] s(3,3)->s(3,3) s(3,3)->s(3,4) s(1,4)->s(1,4) s(2,4)->s(2,4) s(3,4)->s(3,3)
#> [31] s(3,4)->s(2,4) s(3,4)->s(3,4)

plot_transition_graph(Maze)

plot_transition_graph(Maze,
  vertex.size = 20,
  edge.label.cex = .1, edge.arrow.size = .5, margin = .5
)


## Plot using the igraph library
library(igraph)
#> 
#> Attaching package: ‘igraph’
#> The following objects are masked from ‘package:stats’:
#> 
#>     decompose, spectrum
#> The following object is masked from ‘package:base’:
#> 
#>     union
plot(g)


# plot with a different layout
plot(g,
  layout = igraph::layout_with_sugiyama,
  vertex.size = 20,
  edge.label.cex = .6
)


## Use visNetwork (if installed)
if (require(visNetwork)) {
  g_vn <- toVisNetworkData(g)
  nodes <- g_vn$nodes
  edges <- g_vn$edges

  visNetwork(nodes, edges) %>%
    visNodes(physics = FALSE) %>%
    visEdges(smooth = list(type = "curvedCW", roundness = .6), arrows = "to")
}
#> Loading required package: visNetwork

{"x":{"nodes":{"id":["s(1,1)","s(2,1)","s(3,1)","s(1,2)","s(3,2)","s(1,3)","s(2,3)","s(3,3)","s(1,4)","s(2,4)","s(3,4)"],"color":["#FF0000","#FF8B00","#E8FF00","#5DFF00","#00FF2E","#00FFB9","#00B9FF","#002EFF","#5D00FF","#E800FF","#FF008B"],"label":["s(1,1)","s(2,1)","s(3,1)","s(1,2)","s(3,2)","s(1,3)","s(2,3)","s(3,3)","s(1,4)","s(2,4)","s(3,4)"]},"edges":{"from":["s(1,1)","s(1,1)","s(1,1)","s(2,1)","s(2,1)","s(2,1)","s(3,1)","s(3,1)","s(3,1)","s(1,2)","s(1,2)","s(1,2)","s(3,2)","s(3,2)","s(3,2)","s(1,3)","s(1,3)","s(1,3)","s(1,3)","s(2,3)","s(2,3)","s(2,3)","s(2,3)","s(3,3)","s(3,3)","s(3,3)","s(3,3)","s(1,4)","s(2,4)","s(3,4)","s(3,4)","s(3,4)"],"to":["s(1,1)","s(2,1)","s(1,2)","s(1,1)","s(2,1)","s(3,1)","s(2,1)","s(3,1)","s(3,2)","s(1,1)","s(1,2)","s(1,3)","s(3,1)","s(3,2)","s(3,3)","s(1,2)","s(1,3)","s(2,3)","s(1,4)","s(1,3)","s(2,3)","s(3,3)","s(2,4)","s(3,2)","s(2,3)","s(3,3)","s(3,4)","s(1,4)","s(2,4)","s(3,3)","s(2,4)","s(3,4)"],"label":["up (0.9)/\nright (0.1)/\ndown (0.1)/\nleft (0.9)","right (0.1)/\ndown (0.8)/\nleft (0.1)","up (0.1)/\nright (0.8)/\ndown (0.1)","up (0.8)/\nright (0.1)/\nleft (0.1)","up (0.2)/\nright (0.8)/\ndown (0.2)/\nleft (0.8)","right (0.1)/\ndown (0.8)/\nleft (0.1)","up (0.8)/\nright (0.1)/\nleft (0.1)","up (0.1)/\nright (0.1)/\ndown (0.9)/\nleft (0.9)","up (0.1)/\nright (0.8)/\ndown (0.1)","up (0.1)/\ndown (0.1)/\nleft (0.8)","up (0.8)/\nright (0.2)/\ndown (0.8)/\nleft (0.2)","up (0.1)/\nright (0.8)/\ndown (0.1)","up (0.1)/\ndown (0.1)/\nleft (0.8)","up (0.8)/\nright (0.2)/\ndown (0.8)/\nleft (0.2)","up (0.1)/\nright (0.8)/\ndown (0.1)","up (0.1)/\ndown (0.1)/\nleft (0.8)","up (0.8)/\nright (0.1)/\nleft (0.1)","right (0.1)/\ndown (0.8)/\nleft (0.1)","up (0.1)/\nright (0.8)/\ndown (0.1)","up (0.8)/\nright (0.1)/\nleft (0.1)","up (0.1)/\ndown (0.1)/\nleft (0.8)","right (0.1)/\ndown (0.8)/\nleft (0.1)","up (0.1)/\nright (0.8)/\ndown (0.1)","up (0.1)/\ndown (0.1)/\nleft (0.8)","up (0.8)/\nright (0.1)/\nleft (0.1)","right (0.1)/\ndown (0.8)/\nleft (0.1)","up (0.1)/\nright (0.8)/\ndown (0.1)","up/\nright/\ndown/\nleft","up/\nright/\ndown/\nleft","up (0.1)/\ndown (0.1)/\nleft (0.8)","up (0.8)/\nright (0.1)/\nleft (0.1)","up (0.1)/\nright (0.9)/\ndown (0.9)/\nleft (0.1)"]},"nodesToDataframe":true,"edgesToDataframe":true,"options":{"width":"100%","height":"100%","nodes":{"shape":"dot","physics":false},"manipulation":{"enabled":false},"edges":{"arrows":"to","smooth":{"type":"curvedCW","roundness":0.6}}},"groups":null,"width":null,"height":null,"idselection":{"enabled":false},"byselection":{"enabled":false},"main":null,"submain":null,"footer":null,"background":"rgba(0, 0, 0, 0)"},"evals":[],"jsHooks":[]}
```
