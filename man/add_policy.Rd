% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_policy.R
\name{add_policy}
\alias{add_policy}
\title{Add a Policy to an MDP Problem Description}
\usage{
add_policy(model, policy)
}
\arguments{
\item{model}{a \link{MDP} model description.}

\item{policy}{a policy data.frame.}
}
\value{
The model description with the added policy.
}
\description{
Add a policy to an MDP problem description allows the user to
test policies on modified problem descriptions or to test manually created
policies.
}
\details{
The new policy needs to be a data.frame with one row for each state in the
order the states are defined in the model. The only required column is
\itemize{
\item \code{action}: the action prescribed in the state corresponding to the row.
}

Optional columns are
\itemize{
\item \code{state}: the state names in the order of the states in the model.
The needed names can be obtained by from the \verb{$states} element of the model.
\item \code{U}: with the utility given by the value function for the state.
}
}
\examples{
data(Maze)

sol <- solve_MDP(Maze)
sol

policy(sol)
reward(sol)

# Add a random policy
random_pol <- random_policy(Maze)
random_pol
sol_random <- add_policy(Maze, random_pol)
policy(sol_random)
reward(sol_random)
}
\seealso{
Other MDP: 
\code{\link{MDP}()},
\code{\link{accessors}},
\code{\link{act}()},
\code{\link{available_actions}()},
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

Other policy: 
\code{\link{action}()},
\code{\link{bellman_update}()},
\code{\link{policy}()},
\code{\link{policy_evaluation}()},
\code{\link{q_values}()},
\code{\link{reward}()},
\code{\link{value_function}()}
}
\author{
Michael Hahsler
}
\concept{MDP}
\concept{policy}
