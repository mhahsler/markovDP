# Bellman Update and Bellman operator

Update the value function with a Bellman update.

## Usage

``` r
bellman_update(model, V)

bellman_operator(model, pi, V)
```

## Arguments

- model:

  an MDP problem specification.

- V:

  a vector representing the value function. A single 0 can be used as a
  shorthand for a value function with all 0s.

- pi:

  a policy as a data.frame with at least columns for states and action.
  If `NULL`, then the policy in model is used.

## Value

a list with the updated state value vector U and the taken actions pi.

## Details

The Bellman update updates a value function given the model by applying
the Bellman equation as an update rule for each state:

\$\$v\_{k+1}(s) \leftarrow \max\_{a \in \mathcal{A}(s)} \sum\_{s'} p(s'
\| s,a) \[r(s,a, s') + \gamma v_k(s')\]\$\$

The Bellman update moves the estimated value function \\V\\ closer to
the optimal value function \\v\_\*\\.

The Bellman operator \\B\_\pi\\ updates a value function given the
model, and a policy \\\pi\\:

\$\$(B\_\pi v)(s) = \sum\_{a \in \mathcal{A}} \pi(a\|s) \sum\_{s'} p(s'
\| s,a) \[r(s,a,s') + \gamma v(s')\]\$\$

The Bellman error is \\\delta = B\_\pi v - v\\. The Bellman operator
reduces the Bellman error and moves the value function closer to the
fixed point of the true value function:

\$\$v\_\pi = B\_\pi v\_\pi.\$\$

## References

Sutton, R. S., Barto, A. G. (2020). Reinforcement Learning: An
Introduction. Second edition. The MIT Press.

## See also

Other value_function:
[`Q_values()`](http://michael.hahsler.net/markovDP/reference/Q_values.md),
[`value_function()`](http://michael.hahsler.net/markovDP/reference/value_function.md)

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

# single Bellman update from an all-zero value function
bellman_update(Maze, V = 0)
#> $V
#>  [1] -0.04 -0.04 -0.04 -0.04 -0.04  0.76 -0.04 -0.04  0.00  0.00 -0.04
#> 
#> $pi
#> s(1,1) s(2,1) s(3,1) s(1,2) s(3,2) s(1,3) s(2,3) s(3,3) s(1,4) s(2,4) s(3,4) 
#>  right     up  right  right   left  right   left  right     up     up   down 
#> Levels: up right down left
#> 
#> $Q
#>           up right  down  left
#> s(1,1) -0.04 -0.04 -0.04 -0.04
#> s(2,1) -0.04 -0.04 -0.04 -0.04
#> s(3,1) -0.04 -0.04 -0.04 -0.04
#> s(1,2) -0.04 -0.04 -0.04 -0.04
#> s(3,2) -0.04 -0.04 -0.04 -0.04
#> s(1,3)  0.06  0.76  0.06 -0.04
#> s(2,3) -0.14 -0.84 -0.14 -0.04
#> s(3,3) -0.04 -0.04 -0.04 -0.04
#> s(1,4)  0.00  0.00  0.00  0.00
#> s(2,4)  0.00  0.00  0.00  0.00
#> s(3,4) -0.84 -0.14 -0.04 -0.14
#> 

# perform simple value iteration for 10 iterations
V <- 0
for (i in seq(10))
  V <- bellman_update(Maze, V)$V
  
V
#>  [1] 0.8089656 0.7536291 0.6754400 0.8676516 0.5902298 0.9177724 0.6601727
#>  [8] 0.5771594 0.0000000 0.0000000 0.3509586
```
