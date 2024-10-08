% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/available_actions.R
\name{available_actions}
\alias{available_actions}
\title{Available Actions in a State}
\usage{
available_actions(model, state, neg_inf_reward = TRUE, stay_in_place = TRUE)
}
\arguments{
\item{model}{a \link{MDP} object.}

\item{state}{a character vector of length one specifying the state.}

\item{neg_inf_reward}{logical; consider an action that produced \code{-Inf} reward to all
end states unavailable?}

\item{stay_in_place}{logical; consider an action that results in the same state
with a probability of 1 as unavailable. Note that this will mean that
absorbing states have no available action!}
}
\value{
a character vector with the available actions.

a vector with the available actions.
}
\description{
Determine the set of actions available in a state.
}
\details{
Unavailable actions are modeled as actions that have an immediate
reward of \code{-Inf} in the reward function.
}
\examples{
data(DynaMaze)
gw_plot(DynaMaze)

# The the following actions are always available:
DynaMaze$actions

# only right and down is unavailable for s(1,1) because they
#   make the agent stay in place.
available_actions(DynaMaze, state = "s(1,1)")

# An action that leaves the grid currently is allowed but does not do
# anything.
act(DynaMaze, "s(1,1)", "up")
}
\seealso{
Other MDP: 
\code{\link{MDP}()},
\code{\link{accessors}},
\code{\link{act}()},
\code{\link{add_policy}()},
\code{\link{bellman_update}()},
\code{\link{gridworld}},
\code{\link{policy_evaluation}()},
\code{\link{q_values}()},
\code{\link{regret}()},
\code{\link{sample_MDP}()},
\code{\link{solve_MDP}()},
\code{\link{transition_graph}()},
\code{\link{unreachable_and_absorbing}},
\code{\link{value_function}()}
}
\author{
Michael Hahsler
}
\concept{MDP}
