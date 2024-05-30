
# <img src="man/figures/logo.svg" align="right" height="139" /> R package mdp - Infrastructure for Discrete-Time Markov Decision Processes (MDP)

[![r-universe
status](https://mhahsler.r-universe.dev/badges/mdp)](https://mhahsler.r-universe.dev/mdp)
[![Package on
CRAN](http://www.r-pkg.org/badges/version/mdp)](https://CRAN.R-project.org/package=mdp)
[![CRAN RStudio mirror
downloads](http://cranlogs.r-pkg.org/badges/mdp)](https://CRAN.R-project.org/package=mdp)

## Introduction

A Markov decision process (MDP) is a discrete-time stochastic control
process. In each time step, an agent can perform actions which affect
the system (i.e., may cause the system state to change). The agent’s
goal is to maximize its expected future rewards that depend on the
sequence of system state and the agent’s actions in the future. Solving
the MDP means finding the optimal (or at least a good) policy that
guides the agent’s actions.

The package also interfaces to the following algorithms:

- Dynamic Programming (see Russell and Norvig, 2021)

  - **Value Iteration**
  - **Modified Policy Iteration**

- **Linear Programming**

- Termporal Differencing (see Sutton and Barto, 2020)

  - **Q-Learning**
  - **Sarsa**
  - **Expected Sarsa**

To cite package ‘mdp’ in publications use:

> Hahsler M (????). *mdp: Infrastructure for Discrete-Time Markov
> Decision Processes (MDP)*. R package version 0.99.0,
> <https://github.com/mhahsler/mdp>.

    @Manual{,
      title = {mdp: Infrastructure for Discrete-Time Markov Decision Processes (MDP)},
      author = {Michael Hahsler},
      note = {R package version 0.99.0},
      url = {https://github.com/mhahsler/mdp},
    }

## Installation

**Stable CRAN version:** Install from within R with

``` r
install.packages("mdp")
```

**Current development version:** Install from
[r-universe.](https://mhahsler.r-universe.dev/mdp)

``` r
install.packages("mdp",
    repos = c("https://mhahsler.r-universe.dev". "https://cloud.r-project.org/"))
```

## Usage

Solving the simple maze

``` r
library("mdp")
data("Maze")
Maze
```

    ## MDP, list - Stuart Russell's 3x4 Maze
    ##   Discount factor: 1
    ##   Horizon: Inf epochs
    ##   Size: 11 states / 4 actions
    ##   Start: s(3,1)
    ## 
    ##   List components: 'name', 'discount', 'horizon', 'states', 'actions',
    ##     'transition_prob', 'reward', 'info', 'start'

``` r
gridworld_plot(Maze)
```

![](inst/README_files/problem-1.png)<!-- -->

``` r
sol <- solve_MDP(model = Maze)
sol
```

    ## MDP, list - Stuart Russell's 3x4 Maze
    ##   Discount factor: 1
    ##   Horizon: Inf epochs
    ##   Size: 11 states / 4 actions
    ##   Start: s(3,1)
    ##   Solved:
    ##     Method: 'value_iteration'
    ##     Solution converged: TRUE
    ## 
    ##   List components: 'name', 'discount', 'horizon', 'states', 'actions',
    ##     'transition_prob', 'reward', 'info', 'start', 'solution'

Display the value function.

``` r
plot_value_function(sol)
```

![](inst/README_files/value_function-1.png)<!-- -->

``` r
gridworld_plot(sol)
```

![](inst/README_files/gridworld_plot-1.png)<!-- -->

## Acknowledgments

Development of this package was supported in part by National Institute
of Standards and Technology (NIST) under grant number
[60NANB17D180](https://www.nist.gov/ctl/pscr/safe-net-integrated-connected-vehicle-computing-platform).

## References

Russell, S., Norvig, P. (2021). Artificial Intelligence: A Modern
Approach. Fourth edition. Prentice Hall.

Sutton, R. S., Barto, A. G. (2020). Reinforcement Learning: An
Introduction. Second edition. The MIT Press.
