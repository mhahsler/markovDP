#' Linear Function Approximation
#'
#' Approximate a Q-function a value function or a policy using linear
#' function approximation. These approximation functions are 
#' used internally by: [solve_MDP_APPROX]
#'
#' ## Linear Approximation
#' The state-action value function is approximated by
#' \deqn{\hat{q}(s,a) = \boldsymbol{w}^\top\phi(s,a),}
#'
#' where \eqn{\boldsymbol{w} \in \mathbb{R}^n} is a weight vector
#' and \eqn{\phi: S \times A \rightarrow \mathbb{R}^n}  is a
#' feature function that
#' maps each state-action pair to a feature vector.
#' Linear approximation has a single optimum and can be optimized using
#' a simple update rule following the gradient of the state-action function
#' \deqn{\nabla \hat{q}(s,a,\boldsymbol{w}) = \phi(s,a).}
#'
#' Value function approximation works in the same way but \eqn{\phi(s)} is used.
#' 
#' ## State-action Feature Vector Construction
#'
#' For a small number of actions, we can
#' follow the construction described by Geramifard et al (2013)
#' which uses a state feature function \eqn{\phi: S \rightarrow \mathbb{R}^{m}}
#' to construct the complete state-action feature vector.
#' Here, we also add an intercept term.
#' The state-action feature vector has length \eqn{1 + |A| \times m}.
#' It has the intercept and then one component for each action. All these components
#' are set to zero and only the active action component is set to \eqn{\phi(s)},
#' where \eqn{s} is the current state.
#' For example, for the state feature
#' vector \eqn{\phi(s) = (3,4)} and action \eqn{a=2} out of three possible
#' actions \eqn{A = \{1, 2, 3\}}, the complete
#' state-action feature vector is \eqn{\phi(s,a) = (0,0,0,1,3,4,0,0,0)}.
#' Each action component has three entries and the 1 represent the intercept
#' for the state feature vector. The zeros represent the components for the
#' two not chosen actions.
#'
#' The construction of the state-action values is implemented in `add_linear_approx_Q_function()`.
#'
#' The state feature function \eqn{\phi()} starts with raw state feature vector
#' \eqn{\mathbf{x} = (x_1,x_2, ..., x_m)} that
#' are either user-specified or constructed by parsing the state labels of
#' form `s(feature list)`.  Then an optional nonlinear transformation
#' can be performed (see [transformation]).
#'
#' ## Internal Representation
#'
#' All approximations are lists with the elements:
#' * `x(s, a)` ... function to construct from state features, action-state features.
#' * `f(s, a, w)` ... approx. function
#' * `gradient(s, a, w)` ... gradient of f at w
#' * `w` ... the weight vector (initially all 0s)
#' * `transformation` ... a transformation kernel function that is applied to state features in x.
#'
#' ## Prediction
#' `approx_value()` calculates approximate value given the weights in
#'    the model or specified weights.
#'
#' @references
#' Alborz Geramifard, Thomas J. Walsh, Stefanie Tellex, Girish Chowdhary, Nicholas Roy, and Jonathan P. How. 2013. A Tutorial on Linear Function Approximators for Dynamic Programming and Reinforcement Learning. Foundations and Trends in Machine Learning 6(4), December 2013, pp. 375-451. \doi{10.1561/2200000042}
#'
#' Konidaris, G., Osentoski, S., & Thomas, P. 2011. Value Function Approximation in Reinforcement Learning Using the Fourier Basis. Proceedings of the AAAI Conference on Artificial Intelligence, 25(1), 380-385. \doi{10.1609/aaai.v25i1.7903}
#'
#' @examples
#' data(Maze)
#'
#' f_q <- q_approx_linear(Maze)
#' f_q
#' approx_value(f_q, 1, 1, model = Maze)
#'
#' f_v <- v_approx_linear(Maze)
#' f_v
#' approx_value(f_v, 1, model = Maze)
#'
#' # TODO: Implement
#' # f_pi <- pi_approx_linear(Maze)
#' # f_pi
#'
#' @name linear_function_approximation
NULL

#' @rdname linear_function_approximation
#' @param model a [MDP] model with defined state features.
#' @param transformation a transformation function. See [transformation].
#' @param ... further parameters are passed on to the [transformation] function. 
#' @export
q_approx_linear  <- function(model,
                             transformation = transformation_linear_basis,
                             ...) {
  transformation <- transformation(model, ...)
  
  # get state feature format
  state_template <- transformation(start(model, as = "feature")[1L, ])
  n_features <- length(state_template)
  feature_names <- names(state_template) 
  n_A <- length(A(model))
  
  # s are state features. Convert to action-state feature.
  x <- function(s, a) {
    # use a transformation function
    s <- transformation(s)
    # TODO: Remove???
    a <- normalize_action_id(a, model)
    
    # one component per action
    x <- numeric(n_A * n_features)
    
    a_pos <- 1L + (a - 1L) * n_features
    x[a_pos:(a_pos + n_features - 1L)] <- s
    
    x
  }
  
  w_init <- setNames(numeric(n_A * n_features), paste(rep(A(model), each = n_features), feature_names, sep = "."))
  
  structure(
    list(
      f = function(s, a, w)
        sum(w * x(s, a)),
      gradient = function(s, a, w)
        x(s, a),
      x = x,
      transformation = transformation,
      w = w_init
    ),
    class = "q_approx_linear"
  )
  
}

# create a linear approx function that can then be added to a model
#' @rdname linear_function_approximation
#' @export
v_approx_linear  <- function(model,
                             transformation = transformation_linear_basis,
                             ...) {
  transformation <- transformation(model, ...)
  
  # get state feature format
  state_template <- transformation(start(model, as = "feature")[1L, ])
  n_features <- length(state_template)
  feature_names <- names(state_template) 
  
  # s are state features.
  # a is ignored!
  x <- function(s, a = NULL) {
    transformation(s)
  }
  
  w_init <- setNames(numeric(n_features), feature_names)
  
  structure(
    list(
      f = function(s, a = NULL, w)
        sum(w * x(s, a)),
      gradient = function(s, a = NULL, w)
        x(s, a),
      x = x,
      transformation = transformation,
      w = w_init
    ),
    class = "v_approx_linear"
  )
  
}

# create a linear approx function that can then be added to a model
#' @rdname linear_function_approximation
#' @export
pi_approx_linear  <- function(model,
                              transformation = transformation_linear_basis,
                              ...) {
  stop("TODO!")
}


#' @rdname linear_function_approximation
#' @param f a linear approximation object.
#' @param state a state or state features.
#' @param action an action. If `NULL`, then the value for all available actions
#'         will be calculated.
#' @param w a weight vector to be used instead of `f$w`. 
#' @export
approx_value <- function(f,
                         state,
                         action = NULL,
                         # NULL is for v approx.
                         w = NULL,
                         model = NULL) {
  w <- w %||% f$w
  
  if (length(w) != length(f$w))
    stop("w does not have the right number of elements.")
    
  if (!is.matrix(state)) {
    if (is.null(model))
      stop("model is required to convert state labels/ids to state features!")
    state <- normalize_state_features(state, model)
  }
  
  if (!is.numeric(action))
    action <- normalize_action_id(action, model)
    
  f$f(state, action, w)
}