% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reward.R
\name{reward}
\alias{reward}
\alias{reward.MDP}
\title{Calculate the Expected Reward of a Policy}
\usage{
reward(x, ...)

\method{reward}{MDP}(x, start = NULL, epoch = 1L, ...)
}
\arguments{
\item{x}{a solved \link{MDP} object.}

\item{...}{further arguments are passed on.}

\item{start}{specification of the current state (see argument start
in \link{MDP} for details). By default the start state defined in
the model as start is used. Multiple states can be specified as rows in a matrix.}

\item{epoch}{epoch for a finite-horizon solutions.}
}
\value{
\code{reward()} returns a vector of reward values, one for each belief if a matrix is specified.

\item{state}{start state to calculate the reward for. if \code{NULL} then the start
state of model is used.}
}
\description{
This function calculates the expected total reward for a MDP policy
given a start state (distribution). The value is calculated using the value function stored
in the MDP solution.
}
\details{
The reward is typically calculated using the value function
of the solution. If these are not available, then \code{\link[=simulate_MDP]{simulate_MDP()}} is
used instead with a warning.
}
\examples{
data("Maze")
Maze
gridworld_matrix(Maze)

sol <- solve_MDP(Maze)
policy(sol)

# reward for the start state s(3,1) specified in the model
reward(sol)

# reward for starting next to the goal at s(1,3)
reward(sol, start = "s(1,3)")

# expected reward when we start from a random state
reward(sol, start = "uniform")
}
\seealso{
Other policy: 
\code{\link{action}()},
\code{\link{policy}()},
\code{\link{policy_evaluation}()},
\code{\link{q_values}()},
\code{\link{value_function}()}
}
\author{
Michael Hahsler
}
\concept{policy}
