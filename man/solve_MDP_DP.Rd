% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solve_MDP_DP.R
\name{solve_MDP_DP}
\alias{solve_MDP_DP}
\title{Solve MDPs using Dynamic Programming}
\usage{
solve_MDP_DP(
  model,
  method = "VI",
  horizon = NULL,
  discount = NULL,
  n = 1000L,
  error = 0.001,
  k_backups = 10L,
  V = NULL,
  ...,
  matrix = TRUE,
  continue = FALSE,
  verbose = FALSE,
  progress = TRUE
)
}
\arguments{
\item{model}{an MDP problem specification.}

\item{method}{string; one of the following solution methods:
\itemize{
\item \code{'VI'} - value iteration
\item \code{'PI'} - policy iteration
\item \code{'GenPS'}, \code{'PS_error'}, \code{'PS_random'} - prioritized sweeping
}}

\item{horizon}{an integer with the number of epochs for problems with a
finite planning horizon. If set to \code{Inf}, the algorithm continues
running iterations till it converges to the infinite horizon solution. If
\code{NULL}, then the horizon specified in \code{model} will be used.}

\item{discount}{discount factor in range \eqn{(0, 1]}. If \code{NULL}, then the
discount factor specified in \code{model} will be used.}

\item{n}{maximum number of iterations allowed to converge. If the
maximum is reached then the non-converged solution is returned with a
warning.}

\item{error}{value iteration: maximum Bellman error allowed for the
convergence criterion.}

\item{k_backups}{policy iteration: maximum number of Bellman backups used in
the iterative policy evaluation step. Policy evaluation typically converges earlier
with a maximum Bellman error less than \code{error}.}

\item{V}{a vector with initial state values. If
\code{NULL}, then the default of a vector of all 0s (\code{\link[=V_zero]{V_zero()}}) is used.}

\item{...}{further parameters are passed on to the solver function.}

\item{matrix}{logical; if \code{TRUE} then matrices for the transition model and
the reward function are taken from the model first. This can be slow if functions
need to be converted or do not fit into memory if the models are large. If these
components are already matrices, then this is very fast. For \code{FALSE}, the
transition probabilities and the reward is extracted when needed. This is slower,
but removes the time and memory requirements needed to calculate the matrices.}

\item{continue}{logical; show a progress bar with estimated time for completion.}

\item{verbose}{logical or a numeric verbose level; if set to \code{TRUE} or \code{1}, the
function displays the used algorithm parameters and progress information.
Levels \verb{>1} provide more detailed solver output in the R console.}

\item{progress}{logical; show a progress bar with estimated time for completion.}
}
\value{
\code{solve_MDP()} returns an object of class MDP or MDPTF which is a list with the
model specifications (\code{model}), the solution (\code{solution}).
The solution is a list with the elements that depend on the used method. Common
elements are:
\itemize{
\item \code{method} with the name of the used method
\item parameters used.
\item \code{converged} did the algorithm converge (\code{NA}) for finite-horizon problems.
\item \code{policy} a list representing the policy graph. The list only has one
element for converged solutions.
}
}
\description{
Solve MDPs via policy and value iteration.
}
\details{
The following dynamic programming methods
are implemented using the algorithms presented in Russell and Norvig (2010).
\itemize{
\item \strong{(Modified) Policy Iteration} (Howard 1960; Puterman and Shin 1978)
starts with a random policy and iteratively performs
a sequence of
\enumerate{
\item (Approximate) policy evaluation to estimate the value function for the
current policy. Iterative policy evaluation can be approximated by
stopping early after \code{k_backups} iterations
(see \code{\link[=policy_evaluation]{policy_evaluation()}}. In this case the algorithm is called
\emph{modified} policy iteration.
\item Policy improvement is performed by updating the policy to be greedy
(see \code{\link[=greedy_policy]{greedy_policy()}})
with respect to the new value function.
The algorithm stops when it converges to a stable policy (i.e., no changes
between two iterations). Note that the policy typically stabilizes before
the value function converges.
}
\item \strong{Value Iteration} (Bellman 1957) starts with
an arbitrary value function (by default all 0s) and iteratively
updates the value function for each state using the Bellman update
equation (see \code{\link[=bellman_update]{bellman_update()}}).

\deqn{v(s) \leftarrow \max_{a \in \mathcal{A}(s)} \sum_{s'} p(s' | s,a) [r(s,a, s') + \gamma v(s')]}

The iteration
is terminated when the solution converges or the maximum of \code{n} iterations
has been reached.
Approximate convergence is achieved
for discounted problems (with \eqn{\gamma < 1})
when the maximal value function change for any state \eqn{\delta} is
\deqn{\delta \le \frac{error (1-\gamma)}{\gamma}.}
It can be shown that this means
that no state value is more than
\eqn{error} from the value in the optimal value function. For undiscounted
problems, we use \eqn{\delta \le error}.

A greedy policy
is extracted from the final value function. Value iteration can be seen as
policy iteration with policy evaluation truncated to one step.
\item \strong{Prioritized Sweeping} (Moore and Atkeson, 1993; Andre et al., 1997; Li and Littman, 2008)
approximate the optimal value
function by iteratively adjusting one state at a time. While value and policy iteration
sweep in every iteration through all states, prioritized sweeping
updates states in the order given by their priority.
The priority reflects how much a state value may change
given the most recently updated other states that can be directly reached via an action.
This update order often lead to faster convergence compared
to sweeping the whole state space in regular value iteration.

We implement the two priority update strategies described as \strong{PS} and
\strong{GenPS} by Li and Littman (2008).
\itemize{
\item \strong{PS} (Moore and Atkeson, 1993) updates the priority of a state \eqn{H(s)}
using:
\deqn{
       \forall{s \in \mathcal{S}}: H_{t+1}(s)  \leftarrow \begin{cases}
         \max(H_{t}(s), \Delta_t \max_{a \in \mathcal{A}}(p(s_t|s,a)) \text{ for } s \ne s_{t+1} \\
         \Delta_t \max_{a \in A}(p(s_t|s,a) \text{ for } s = s_{t+1}
         \end{cases}
     }

where \eqn{\Delta_t = |V_{t+1}(s_t) - V_t(s_t)| = |E(s_t; V_{t+1})|}, i.e.,
the Bellman error for the updated state.
\item \strong{GenPS} (Andre et al., 1997) updates all state priorities using their
current Bellman error:

\deqn{\forall{s \in \mathcal{S}}: H_{t+1}(s) \leftarrow |E(s; V_{t+1})|}

where \eqn{E(s; V_{t+1}) = \max_{a \in A} \left[R(s,a) + \gamma \sum_{s \in S} p(s'|s,a) V(s')\right] - V(s)}
is a state's Bellman error.
}

The update method can be chosen using the additional parameter \code{H_update}
as the character string \code{"PS_random"}, \code{"PS_error"} or \code{"GenPS"}.
The default is \code{H_update = "GenPS"}. For PS, random means that the
priority vector is initialized with random values (larger than 0),
and error means they are initialized with the Bellman error as in
GenPS. However, this requires one complete sweep over all states.

This implementation stops updating when the largest priority values
over all states is less than the specified \code{error}.

Since the algorithm does not sweep through the whole state space for each
iteration, \code{n} is converted into an equivalent number of state updates
\eqn{n = n\ |S|}.
}
}
\examples{
data(Maze)

maze_solved <- solve_MDP(Maze, method = "DP:VI", verbose = TRUE)
policy(maze_solved)

# use prioritized sweeping (which is known to be fast for mazes)
maze_solved <- solve_MDP(Maze, method = "DP:GenPS", verbose = TRUE)
policy(maze_solved)

# finite horizon
maze_solved <- solve_MDP(Maze, method = "DP:VI", horizon = 3)
policy(maze_solved)
gw_plot(maze_solved, epoch = 1)
gw_plot(maze_solved, epoch = 2)
gw_plot(maze_solved, epoch = 3)

}
\references{
Andre, D., Friedman, N., and Parr, R. 1997. "Generalized prioritized sweeping." In Advances in Neural Information Processing Systems 10, pp. 1001-1007. \href{https://proceedings.neurips.cc/paper_files/paper/1997/file/7b5b23f4aadf9513306bcd59afb6e4c9-Paper.pdf}{NeurIPS Proceedings}

Bellman, Richard. 1957. "A Markovian Decision Process." Indiana University Mathematics Journal 6: 679-84. \url{https://www.jstor.org/stable/24900506}.

Howard, R. A. 1960. "Dynamic Programming and Markov Processes." Cambridge, MA: MIT Press.

Li, Lihong, and Michael Littman. 2008. "Prioritized Sweeping Converges to the Optimal Value Function." DCS-TR-631. Rutgers University. \doi{10.7282/T3TX3JSX}

Moore, Andrew, and C. G. Atkeson. 1993. "Prioritized Sweeping: Reinforcement Learning with Less Data and Less Real Time." Machine Learning 13 (1): 103–30. \doi{10.1007/BF00993104}.

Puterman, Martin L., and Moon Chirl Shin. 1978. "Modified Policy Iteration Algorithms for Discounted Markov Decision Problems." Management Science 24: 1127-37. \doi{10.1287/mnsc.24.11.1127}.

Russell, Stuart J., and Peter Norvig. 2020. Artificial Intelligence: A Modern Approach (4th Edition). Pearson. \url{http://aima.cs.berkeley.edu/}.
}
\seealso{
Other solver: 
\code{\link{schedule}},
\code{\link{solve_MDP}()},
\code{\link{solve_MDP_APPROX}()},
\code{\link{solve_MDP_LP}()},
\code{\link{solve_MDP_MC}()},
\code{\link{solve_MDP_PG}()},
\code{\link{solve_MDP_SAMP}()},
\code{\link{solve_MDP_TD}()}
}
\concept{solver}
