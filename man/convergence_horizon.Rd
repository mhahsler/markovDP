% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convergence_horizon.R
\name{convergence_horizon}
\alias{convergence_horizon}
\title{Estimate the Convergence Horizon for an Infinite-Horizon MDP}
\usage{
convergence_horizon(model, delta = 0.001, n_updates = 10)
}
\arguments{
\item{model}{an MDP model.}

\item{delta}{maximum update error.}

\item{n_updates}{integer; average number of time each state is updated.}
}
\value{
An estimated convergence horizon.
}
\description{
Many sampling-based methods require a finite horizon. Fo infinite horizons,
discounting leads to convergences during a finite horizon. This function
estimates the number of steps till convergence using rules of thumb.
}
\details{
The horizon is estimated differently for the discounted and the undiscounted
case.
\subsection{Discounted Case}{

The effect of the largest reward \eqn{R_{\mathrm{max}}}
update decreases with \eqn{t} as \eqn{\delta_t = \gamma^t R_{\mathrm{max}}}.
The convergence horizon is estimated as the smallest \eqn{t} for which
\eqn{\delta_t < \delta}.
}

\subsection{Undiscounted Case}{

For the undiscounted case, episodes end when an absorbing state is reached.
It cannot be guaranteed that a model will reach an absorbing state.
To avoid infinite loops, we set the maximum horizon such that each entry in
the Q-table is on average updated \code{n_updates} times. This is a very rough
rule ot thumb.
}
}
\examples{
data(Maze)
Maze

convergence_horizon(Maze)

# make the Maze into a discounted problem where future rewards count less.
Maze_discounted <- Maze
Maze_discounted$discount <- .9
Maze_discounted

convergence_horizon(Maze_discounted)
}
\seealso{
Other MDP: 
\code{\link{MDP}()},
\code{\link{Q_values}()},
\code{\link{absorbing_states}()},
\code{\link{act}()},
\code{\link{available_actions}()},
\code{\link{bellman_update}()},
\code{\link{greedy_action}()},
\code{\link{gridworld}},
\code{\link{policy_evaluation}()},
\code{\link{reachable_states}()},
\code{\link{regret}()},
\code{\link{sample_MDP}()},
\code{\link{solve_MDP}()},
\code{\link{start}()},
\code{\link{transition_graph}()},
\code{\link{transition_matrix}()},
\code{\link{unreachable_states}()},
\code{\link{value_function}()}
}
\author{
Michael Hahsler
}
\concept{MDP}
