
# <img src="man/figures/logo.svg" align="right" height="139" /> R package markovDP - Infrastructure for Discrete-Time Markov Decision Processes (MDP)

[![r-universe
status](https://mhahsler.r-universe.dev/badges/markovDP)](https://mhahsler.r-universe.dev/markovDP)
[![Package on
CRAN](http://www.r-pkg.org/badges/version/markovDP)](https://CRAN.R-project.org/package=markovDP)
[![CRAN RStudio mirror
downloads](http://cranlogs.r-pkg.org/badges/markovDP)](https://CRAN.R-project.org/package=markovDP)

## Introduction

A Markov decision process (MDP) ([Bellman 1957](#ref-Bellman1957);
[Howard 1960](#ref-Howard1960)) is a discrete-time stochastic control
process. In each time step, an agent can perform actions which affect
the system (i.e., may cause the system state to change). The agent’s
goal is to maximize its expected future rewards that depend on the
sequence of system state and the agent’s actions in the future. Solving
the MDP means finding the optimal (or at least a good) policy that
guides the agent’s actions.

The `markovDP` package provides the infrastructure to work with MDPs in
R. The focus is on convenience in formulating MDPs in multiple ways, the
support of sparse representations (using sparse matrices, lists and
data.frames) and visualization of results. Some key components are
implemented in C++ to speed up computation. It also provides to the
following popular solving procedures:

- **Dynamic Programming**
  - Value Iteration ([Bellman 1957](#ref-Bellman1957))
  - Modified Policy Iteration ([Howard 1960](#ref-Howard1960); [Puterman
    and Shin 1978](#ref-Puterman1978))
  - Prioritized Sweeping ([Moore and Atkeson 1993](#ref-Moore1993))
- **Linear Programming**
  - Primal Formulation ([Manne 1960](#ref-Manne1960))
- **Termporal Differencing**
  - Q-Learning ([Watkins and Dayan 1992](#ref-Watkins1992))
  - Sarsa ([Rummery and Niranjan 1994](#ref-Rummery1994))
  - Expected Sarsa ([Sutton and Barto 2018](#ref-Sutton1998))

These implementations follow the description is ([Russell and Norvig
2020](#ref-Russell2020)) and ([Sutton and Barto 2018](#ref-Sutton1998)).

To cite package ‘markovDP’ in publications use:

> Hahsler M (????). *markovDP: Infrastructure for Discrete-Time Markov
> Decision Processes (MDP)*. R package version 0.99.0,
> <https://github.com/mhahsler/markovDP>.

    @Manual{,
      title = {markovDP: Infrastructure for Discrete-Time Markov Decision Processes (MDP)},
      author = {Michael Hahsler},
      note = {R package version 0.99.0},
      url = {https://github.com/mhahsler/markovDP},
    }

## Installation

**Current development version:** Install from
[r-universe.](https://mhahsler.r-universe.dev/markovDP)

``` r
install.packages("markovDP",
    repos = c("https://mhahsler.r-universe.dev". "https://cloud.r-project.org/"))
```

## Usage

Solving the simple maze from ([Russell and Norvig
2020](#ref-Russell2020)).

``` r
library("markovDP")
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

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Bellman1957" class="csl-entry">

Bellman, Richard. 1957. “A Markovian Decision Process.” *Indiana
University Mathematics Journal* 6: 679–84.
<https://www.jstor.org/stable/24900506>.

</div>

<div id="ref-Howard1960" class="csl-entry">

Howard, R. A. 1960. *Dynamic Programming and Markov Processes*.
Cambridge, MA: MIT Press.

</div>

<div id="ref-Manne1960" class="csl-entry">

Manne, Alan. 1960. “On the Job-Shop Scheduling Problem.” *Operations
Research* 8 (2): 219–23. <https://doi.org/10.1287/opre.8.2.219>.

</div>

<div id="ref-Moore1993" class="csl-entry">

Moore, Andrew, and C. G. Atkeson. 1993. “Prioritized Sweeping:
Reinforcement Learning with Less Data and Less Real Time.” *Machine
Learning* 13 (1): 103–30. <https://doi.org/10.1007/BF00993104>.

</div>

<div id="ref-Puterman1978" class="csl-entry">

Puterman, Martin L., and Moon Chirl Shin. 1978. “Modified Policy
Iteration Algorithms for Discounted Markov Decision Problems.”
*Management Science* 24: 1127–37.
<https://doi.org/10.1287/mnsc.24.11.1127>.

</div>

<div id="ref-Rummery1994" class="csl-entry">

Rummery, G., and Mahesan Niranjan. 1994. “On-Line Q-Learning Using
Connectionist Systems.” Techreport CUED/F-INFENG/TR 166. Cambridge
University Engineering Department.

</div>

<div id="ref-Russell2020" class="csl-entry">

Russell, Stuart J., and Peter Norvig. 2020. *Artificial Intelligence: A
Modern Approach (4th Edition)*. Pearson. <http://aima.cs.berkeley.edu/>.

</div>

<div id="ref-Sutton1998" class="csl-entry">

Sutton, Richard S., and Andrew G. Barto. 2018. *Reinforcement Learning:
An Introduction*. Second. The MIT Press.
<http://incompleteideas.net/book/the-book-2nd.html>.

</div>

<div id="ref-Watkins1992" class="csl-entry">

Watkins, Christopher J. C. H., and Peter Dayan. 1992. “Q-Learning.”
*Machine Learning* 8 (3): 279–92. <https://doi.org/10.1007/BF00992698>.

</div>

</div>
