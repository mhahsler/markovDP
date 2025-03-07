% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solve_MDP_APPROX.R
\name{solve_MDP_APPROX}
\alias{solve_MDP_APPROX}
\alias{add_linear_approx_Q_function}
\alias{add_linear_approx_Q_function.MDP}
\alias{add_linear_approx_Q_function.MDPTF}
\alias{create_basis_coefs}
\alias{transformation_linear_basis}
\alias{transformation_polynomial_basis}
\alias{transformation_RBF_basis}
\alias{transformation_fourier_basis}
\alias{approx_Q_value}
\alias{approx_greedy_action}
\alias{approx_greedy_policy}
\alias{approx_V_plot}
\title{Episodic Semi-gradient Sarsa with Linear Function Approximation}
\usage{
solve_MDP_APPROX(
  model,
  method = "semi_gradient_sarsa",
  horizon = NULL,
  discount = NULL,
  alpha = schedule_exp(0.2, 0.1),
  epsilon = schedule_exp(1, 0.1),
  n,
  w = NULL,
  ...,
  matrix = TRUE,
  continue = FALSE,
  progress = TRUE,
  verbose = FALSE
)

add_linear_approx_Q_function(model, ...)

\method{add_linear_approx_Q_function}{MDP}(
  model,
  state_features = NULL,
  transformation = transformation_linear_basis,
  ...
)

\method{add_linear_approx_Q_function}{MDPTF}(
  model,
  state_features = NULL,
  transformation = transformation_linear_basis,
  ...
)

create_basis_coefs(dim, order)

transformation_linear_basis(model, min = NULL, max = NULL, intercept = TRUE)

transformation_polynomial_basis(
  model,
  min = NULL,
  max = NULL,
  order,
  coefs = NULL
)

transformation_RBF_basis(
  model,
  min = NULL,
  max = NULL,
  n,
  centers = NULL,
  var = NULL
)

transformation_fourier_basis(
  model,
  min = NULL,
  max = NULL,
  order,
  coefs = NULL
)

approx_Q_value(model, state = NULL, action = NULL, w = NULL)

approx_greedy_action(model, state, w = NULL, epsilon = 0)

approx_greedy_policy(model, w = NULL)

approx_V_plot(
  model,
  min = NULL,
  max = NULL,
  w = NULL,
  res = 25,
  col = hcl.colors(res, "YlOrRd", rev = TRUE),
  image = TRUE,
  contour = TRUE,
  main = NULL,
  ...
)
}
\arguments{
\item{model}{an MDP problem specification.}

\item{method}{string; one of the following solution methods: \code{'semi_gradient_sarsa'}}

\item{horizon}{an integer with the number of epochs for problems with a
finite planning horizon. If set to \code{Inf}, the algorithm continues
running iterations till it converges to the infinite horizon solution. If
\code{NULL}, then the horizon specified in \code{model} will be used.}

\item{discount}{discount factor in range \eqn{(0, 1]}. If \code{NULL}, then the
discount factor specified in \code{model} will be used.}

\item{alpha}{step size (learning rate). A scalar value between 0 and 1 or a
\link{schedule}.}

\item{epsilon}{used for the \eqn{\epsilon}-greedy behavior policies. A scalar value between 0 and 1 or a
\link{schedule}.}

\item{n}{number of episodes used for learning.}

\item{w}{a weight vector}

\item{...}{further parameters are passed on to the solver function.}

\item{matrix}{logical; if \code{TRUE} then matrices for the transition model and
the reward function are taken from the model first. This can be slow if functions
need to be converted or do not fit into memory if the models are large. If these
components are already matrices, then this is very fast. For \code{FALSE}, the
transition probabilities and the reward is extracted when needed. This is slower,
but removes the time and memory requirements needed to calculate the matrices.}

\item{continue}{logical; show a progress bar with estimated time for completion.}

\item{progress}{logical; show a progress bar with estimated time for completion.}

\item{verbose}{logical or a numeric verbose level; if set to \code{TRUE} or \code{1}, the
function displays the used algorithm parameters and progress information.
Levels \verb{>1} provide more detailed solver output in the R console.}

\item{state_features}{a matrix with state features. Each row is the feature
vector for a state.}

\item{transformation}{a transformation function that is applied to the feature vector
before it is used in the linear approximation.}

\item{dim}{number of features to describe a state.}

\item{order}{order for the Fourier basis.}

\item{min, max}{vectors with the minimum and maximum values for each feature.
This is used to scale the feature to the \eqn{[0,1]} interval for the
Fourier basis.}

\item{intercept}{logical; add an intercept term to the linear basis?}

\item{coefs}{an optional matrix or data frame to specify the set of
coefficient values for the
Fourier basis (overrides \code{order}).}

\item{centers}{a matrix with the centers for the RBF. By default a regular
grid with n steps per feature dimension is used.}

\item{var}{a scalar with the variance used for the RBF.}

\item{state}{a state (index or name)}

\item{action}{an action (index or name)}

\item{res}{resolution as the number of values sampled from each feature.}

\item{col}{colors for the passed on to \code{\link[=image]{image()}}.}

\item{image, contour}{logical; include the false color image or the
contours in the plot?}

\item{main}{title for the plot.}
}
\value{
\code{solve_MDP()} returns an object of class MDP or MDPTF which is a list with the
model specifications (\code{model}), the solution (\code{solution}).
The solution is a list with the elements:
\itemize{
\item \code{policy} a list representing the policy graph. The list only has one
element for converged solutions.
\item \code{converged} did the algorithm converge (\code{NA}) for finite-horizon problems.
\item \code{delta} final \eqn{\delta} (value iteration and infinite-horizon only)
\item \code{iterations} number of iterations to convergence (infinite-horizon only)
}
}
\description{
Solve the MDP control problem using state-value approximation
by semi-gradient Sarsa (temporal differencing) for
episodic problems.
}
\details{
\subsection{Linear Approximation}{

The state-action value function is approximated by
\deqn{\hat{q}(s,a) = \boldsymbol{w}^\top\phi(s,a),}

where \eqn{\boldsymbol{w} \in \mathbb{R}^n} is a weight vector
and \eqn{\phi: S \times A \rightarrow \mathbb{R}^n}  is a
feature function that
maps each state-action pair to a feature vector.
Linear approximation has a single optimum and can be optimized using
a simple update rule following the gradient of the state-action function
\deqn{\nabla \hat{q}(s,a,\boldsymbol{w}) = \phi(s,a).}
}

\subsection{State-action Feature Vector Construction}{

For a small number of actions, we can
follow the construction described by Geramifard et al (2013)
which uses a state feature function \eqn{\phi: S \rightarrow \mathbb{R}^{m}}
to construct the complete state-action feature vector.
Here, we also add an intercept term.
The state-action feature vector has length \eqn{1 + |A| \times m}.
It has the intercept and then one component for each action. All these components
are set to zero and only the active action component is set to \eqn{\phi(s)},
where \eqn{s} is the current state.
For example, for the state feature
vector \eqn{\phi(s) = (3,4)} and action \eqn{a=2} out of three possible
actions \eqn{A = \{1, 2, 3\}}, the complete
state-action feature vector is \eqn{\phi(s,a) = (0,0,0,1,3,4,0,0,0)}.
Each action component has three entries and the 1 represent the intercept
for the state feature vector. The zeros represent the components for the
two not chosen actions.

The construction of the state-action values is implemented in \code{add_linear_approx_Q_function()}.

The state feature function \eqn{\phi()} starts with raw state feature vector
\eqn{\mathbf{x} = (x_1,x_2, ..., x_m)} that
are either user-specified or constructed by parsing the state labels of
form \verb{s(feature list)}.  Then an optional transformations called basis functions
can be applied. Implemented basis functions are:
\itemize{
\item Linear: no additional transformation is applied giving \eqn{\phi_0(s) = 1}
for the intercept and \eqn{\phi_i(s) = x_i} for \eqn{i = \{1, 2, ..., m\}}.
\item Polynomial basis: \eqn{\phi_i(s) = \prod_{j=1}^m x_j^{c_{i,j}}}, where
\eqn{c_{1,j}} is an integer between 0 and \eqn{n} for and order \eqn{n} polynomial basis.
\item Fourier basis: \eqn{\phi_i(s) = \text{cos}(\pi\mathbf{c}^i \cdot \mathbf{x})},
where \eqn{\mathbf{c}^i = [c_1, c_2, ..., c_m]} with \eqn{c_j = [0, ..., n]}, where
\eqn{n} is the order of the basis. The components of the feature vector
\eqn{x} are assumed to be scaled to the interval \eqn{[0,1]}. The fourier
basis transformation is implemented in \code{transformation_fourier_basis()}.
\code{min} and \code{max} are the minimums and maximums for each feature vector component
used to resale them to \eqn{[0,1]} using \eqn{\frac{x_i - min_i}{max_i - min_i}}

Details of this transformation are described in Konidaris et al (2011).
}
}

\subsection{Helper Functions}{

The following helper functions for using approximation are available:
\itemize{
\item \code{approx_Q_value()} calculates approximate Q values given the weights in
the model or specified weights.
\item \code{approx_greedy_action()} uses approximate Q values given the weights in
the model or specified weights to find the the greedy action for a
state.
\item \code{approx_greedy_policy()} calculates the greedy-policy
for the approximate Q values given the weights in
the model or specified weights.
}
}

\subsection{Episodic Semi-gradient Sarsa}{

The implementation follows the temporal difference algorithm
episodic Semi-gradient Sarsa algorithm given in Sutton and Barto (2018).
}

\subsection{Schedules}{
\itemize{
\item epsilon schedule: \code{t} is increased by each processed episode.
\item alpha schedule: \code{t} is increased by each processed episode.
}
}
}
\examples{
# Example 1: A maze without walls. The step cost is 1. The start is top-left and
# the goal (+100 reward) is bottom-right.
# This is the ideal problem for a linear approximation of the Q-function
# using the x/y location as state features.

m <- gw_maze_MDP(c(5, 5), start = "s(1,1)", goal = "s(5,5)")

# construct state features as the x/y coordinates in the gridworld
state_features <- gw_s2rc(S(m))
state_features

m <- add_linear_approx_Q_function(m, state_features)

# constructed state-action features (X) and approximate Q function
# and gradient
m$approx_Q_function

set.seed(1000)
sol <- solve_MDP_APPROX(m, horizon = 1000, n = 100)

gw_plot(sol)
gw_matrix(sol, what = "value")

# learned weights and state values
sol$solution$w

# the approximate value function can be visualized for states 
# with two features.
approx_V_plot(sol)

# extracting approximate Q-values
approx_greedy_action(sol, "s(4,5)")
approx_Q_value(sol, "s(4,5)", "down")
approx_Q_value(sol)

# extracting a greedy policy using the approximate Q-values
approx_greedy_policy(sol)


# Example 2: Stuart Russell's 3x4 Maze using linear basis approximation
# The wall and the -1 absorbing state make linear approximation
# using just the position directly more difficult.

data(Maze)
gw_plot(Maze)

# if no state features are specified, then they are constructed
# by parsing the state label of the form s(feature list).
Maze_approx <- add_linear_approx_Q_function(Maze)

set.seed(1000)
sol <- solve_MDP_APPROX(Maze_approx, horizon = 100, n = 100,
                     alpha = schedule_harmonic(start = .2, n = n), 
                     epsilon = schedule_harmonic(start = 10, n = n))
gw_plot(sol)
gw_matrix(sol, what = "value")
approx_V_plot(sol, res = 20)


# Example 3: Stuart Russell's 3x4 Maze using
#            order-2 Fourier basis for approximation

Maze_approx <- add_linear_approx_Q_function(Maze,
      transformation = transformation_fourier_basis, order = 2)

set.seed(1000)
sol <- solve_MDP_APPROX(Maze_approx, horizon = 100, n = 100,
                    alpha = schedule_exp(0.2, .1),
                    epsilon = schedule_exp(1, .1))
                    
gw_plot(sol)
gw_matrix(sol, what = "value")
approx_V_plot(sol, res = 20)


# Example 4: Stuart Russell's 3x4 Maze using
#            order-2 polynomial basis for approximation

Maze_approx <- add_linear_approx_Q_function(Maze,
      transformation = transformation_polynomial_basis, order = 2)

set.seed(1000)
sol <- solve_MDP_APPROX(Maze_approx, horizon = 100, n = 100,
                    alpha = schedule_exp(0.2, .01),
                    epsilon = schedule_exp(1, .1))
  
gw_plot(sol)
gw_matrix(sol, what = "value")
approx_V_plot(sol, res = 20)

# Example 5: Stuart Russell's 3x4 Maze using RBF with
#            4^2 centers for approximation

Maze_approx <- add_linear_approx_Q_function(Maze,
      transformation = transformation_RBF_basis, n = 4)

set.seed(1000)
sol <- solve_MDP_APPROX(Maze_approx, horizon = 100, n = 100,
                    alpha = schedule_exp(0.2, .001),
                    epsilon = schedule_exp(1, .1))
gw_plot(sol)
gw_matrix(sol, what = "value")
approx_V_plot(sol, res = 20)

}
\references{
Sutton, Richard S., and Andrew G. Barto. 2018. Reinforcement Learning: An Introduction. Second. The MIT Press. \url{http://incompleteideas.net/book/the-book-2nd.html}.

Alborz Geramifard, Thomas J. Walsh, Stefanie Tellex, Girish Chowdhary, Nicholas Roy, and Jonathan P. How. 2013. A Tutorial on Linear Function Approximators for Dynamic Programming and Reinforcement Learning. Foundations and Trends in Machine Learning 6(4), December 2013, pp. 375-451. \doi{10.1561/2200000042}

Konidaris, G., Osentoski, S., & Thomas, P. 2011. Value Function Approximation in Reinforcement Learning Using the Fourier Basis. Proceedings of the AAAI Conference on Artificial Intelligence, 25(1), 380-385. \doi{10.1609/aaai.v25i1.7903}
}
\seealso{
Other solver: 
\code{\link{schedule}},
\code{\link{solve_MDP}()},
\code{\link{solve_MDP_APPROX_LAMBDA}()},
\code{\link{solve_MDP_DP}()},
\code{\link{solve_MDP_LP}()},
\code{\link{solve_MDP_SAMP}()},
\code{\link{solve_MDP_TD}()}

Other MDPTF: 
\code{\link{MDPTF}()},
\code{\link{absorbing_states}()},
\code{\link{act}()},
\code{\link{reachable_states}()},
\code{\link{sample_MDP.MDPTF}()},
\code{\link{solve_MDP}()},
\code{\link{solve_MDP_APPROX_LAMBDA}()},
\code{\link{start}()}
}
\concept{MDPTF}
\concept{solver}
