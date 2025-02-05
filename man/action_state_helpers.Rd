% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/action_state_helpers.R
\name{action_state_helpers}
\alias{action_state_helpers}
\alias{normalize_state}
\alias{normalize_state_id}
\alias{normalize_state_label}
\alias{normalize_action}
\alias{normalize_action_id}
\alias{normalize_action_label}
\alias{state2features}
\alias{features2state}
\alias{s}
\alias{normalize_state_features}
\title{Conversions for Action and State IDs and Labels}
\usage{
normalize_state(state, model)

normalize_state_id(state, model)

normalize_state_label(state, model)

normalize_action(action, model)

normalize_action_id(action, model)

normalize_action_label(action, model)

state2features(state)

features2state(x)

s(...)

normalize_state_features(state, model = NULL)
}
\arguments{
\item{state}{a state labels}

\item{model}{a MDP model}

\item{action}{a action labels}

\item{x}{a state feature vector or a matrix of state feature vectors as rows.}

\item{...}{features that should be converted into a row
vector used to describe a state.}
}
\value{
Functions ending in
\itemize{
\item \verb{_id} return an integer id,
\item \verb{_label} return a character string,
\item \verb{_features} return a state feature matrix,
\item no ending return a factor.
}

Other functions:
\itemize{
\item \code{state2features()} returns a feature vector/matrix.
\item \code{features2state(x)} returns a state label in the format \verb{s(feature list)}.
\item \code{s()} returns a state features row vector.
}
}
\description{
Several helper functions to convert state and action (integer) IDs
to labels and vice versa.
}
\details{
\code{normalize_state()} and \code{normalize_action()} convert labels or ids into a
factor which is the standard representation. If only the label or the
integer id (i.e., the index) is needed, the additional functions can be used.
These are typically a lot faster.

To support a factored state representation as feature vectors,
\code{state2features()} and \code{feature2states()} are provided.
A factored state is represented as a \strong{row} vector for a single
state (ceonveniently created via \code{s()}) or a matrix with row vectors for a
set of states are used. State labels
are constructed in the form \code{s(feature1, feature2, ...)}.
Factored state representation
is used for value function approximation (see \code{\link[=solve_MDP_APPROX]{solve_MDP_APPROX()}}) and
for \link{MDPE} to describe MDP's via a transition function between factored
states.
}
