#' Transformation Functions for Linear Function Approximation
#'
#' Several popular transformation functions applied to 
#' state features used in linear function approximation for 
#' [solve_MDP_APPROX()].
#'
#' The state feature function \eqn{\phi()} uses
#' the raw state feature vectors
#' \eqn{\mathbf{x} = (x_1,x_2, ..., x_m)}
#' which is either user-specified or constructed by parsing the state labels of
#' form `s(feature list)` and then applies    
#' a transformation functions called basis functions.
#' Implemented basis functions are:
#'
#' * Linear: no additional transformation is applied giving \eqn{\phi_0(s) = 1}
#'   for the intercept and \eqn{\phi_i(s) = x_i} for \eqn{i = \{1, 2, ..., m\}}.
#'
#' * Polynomial basis: 
#'   \deqn{\phi_i(s) = \prod_{j=1}^m x_j^{c_{i,j}},} 
#'   where \eqn{c_{1,j}} is an integer between 0 and \eqn{n} 
#'   for and order \eqn{n} polynomial basis.
#'   
#' * Radial Basis: RBF.
#'
#' * Fourier basis: 
#'   \deqn{\phi_i(s) = \text{cos}(\pi\mathbf{c}^i \cdot \mathbf{x}),}
#'   where \eqn{\mathbf{c}^i = [c_1, c_2, ..., c_m]} with \eqn{c_j = [0, ..., n]}, where
#'   \eqn{n} is the order of the basis. The components of the feature vector
#'   \eqn{x} are assumed to be scaled to the interval \eqn{[0,1]}. The fourier
#'   basis transformation is implemented in `transformation_fourier_basis()`.
#'   `min` and `max` are the minimums and maximums for each feature vector component
#'   used to resale them to \eqn{[0,1]} using \eqn{\frac{x_i - min_i}{max_i - min_i}}
#'
#'   Details of this transformation are described in Konidaris et al (2011).
#'
#' @name transformation
#' @aliases transformation
#' 
#' @references
#' Sutton, Richard S., and Andrew G. Barto. 2018. Reinforcement Learning: An Introduction. Second. The MIT Press. [http://incompleteideas.net/book/the-book-2nd.html](http://incompleteideas.net/book/the-book-2nd.html).
#'
#' Alborz Geramifard, Thomas J. Walsh, Stefanie Tellex, Girish Chowdhary, Nicholas Roy, and Jonathan P. How. 2013. A Tutorial on Linear Function Approximators for Dynamic Programming and Reinforcement Learning. Foundations and Trends in Machine Learning 6(4), December 2013, pp. 375-451. \doi{10.1561/2200000042}
#'
#' Konidaris, G., Osentoski, S., & Thomas, P. 2011. Value Function Approximation in Reinforcement Learning Using the Fourier Basis. Proceedings of the AAAI Conference on Artificial Intelligence, 25(1), 380-385. \doi{10.1609/aaai.v25i1.7903}
#'
#'
#' @returns A transformation function
#'
#' @seealso [solve_MDP_APPROX()]
#'
#' @param model the [MDP] model.
#' @param min,max vectors with the minimum and maximum values for each feature.
#'    This is used to scale the feature to the \eqn{[0,1]} interval for the
#'    Fourier basis.
#' @param intercept logical; add an intercept term to the linear basis?
#' @export
transformation_linear_basis <- function(model,
                                        min = NULL,
                                        max = NULL,
                                        intercept = TRUE) {
  rng <- get_state_feature_range(model, min, max)
  min <- rng[1, ]
  max <- rng[2, ]
  
  function(x) {
    # linear basis
    x <- (x - min) / (max - min)
    
    if (intercept)
      x <- c(x0 = 1, x)
    x
  }
}


#' @rdname transformation
#' @export
transformation_polynomial_basis <- function(model,
                                            min = NULL,
                                            max = NULL,
                                            order,
                                            coefs = NULL) {
  rng <- get_state_feature_range(model, min, max)
  min <- rng[1, ]
  max <- rng[2, ]
  dim_s <- length(max)
  
  if (is.null(coefs))
    coefs <- create_basis_coefs(dim_s, order)
  
  function(x) {
    # polynomial basis 
    x <- (x - min) / (max - min)
    
    apply(
      coefs,
      MARGIN = 1,
      FUN = function(c)
        prod(x^c)
    )
  }
}

#' @rdname transformation
#' @param centers a matrix with the centers for the RBF. By default a regular 
#'               grid with n steps (see below) per feature dimension is used.
#' @param var a scalar with the variance used for the RBF.
#' @param n number of centers.
#' @export
transformation_RBF_basis <- function(model,
                                     min = NULL,
                                     max = NULL,
                                     n,
                                     centers = NULL,
                                     var = NULL) {
  rng <- get_state_feature_range(model, min, max)
  min <- rng[1, ]
  max <- rng[2, ]
  dim_s <- length(max)
  
  # create n^d centers
  if (!is.null(centers)) {
    n <- nrow(centers)
  } else {
    if (is.null(n))
      stop("either n or centers has to be specified")
    # create equally spaced centers
    step_size <- 1 / n
    centers <- expand.grid(replicate(dim_s, seq(
      step_size / 2, 1 - (step_size / 2), length.out = n
    ) , simplify = FALSE))
  }
  
  var <- var %||% 2 / (n - 1)
  
  function(x) {
    # RBF basis 
    x <- (x - min) / (max - min)
    apply(
      centers,
      MARGIN = 1,
      FUN = function(c)
        1 / (2 * pi * var)^.5 * exp(-1 * sum((c - x)^2) / 2 / var)
    )
  }
}

#' @rdname transformation
#' @param order order for the Fourier basis.
#' @param coefs an optional matrix or data frame to specify the set of
#'    coefficient values for the
#'    Fourier basis (overrides `order`).
#' @export
transformation_fourier_basis <- function(model,
                                         min = NULL,
                                         max = NULL,
                                         order,
                                         coefs = NULL) {
  rng <- get_state_feature_range(model, min, max)
  min <- rng[1, ]
  max <- rng[2, ]
  dim_s <- length(max)
  
  coefs <- coefs %||% create_basis_coefs(dim_s, order)
  
  function(x) {
    # Fourier basis
    x <- (x - min) / (max - min)
    
    apply(
      coefs,
      MARGIN = 1,
      FUN = function(c)
        cos(pi * x %*% c)
    )
  }
}


#' @rdname transformation
#' @param dim number of features to describe a state.
#' @export
create_basis_coefs <- function(dim, order) {
  as.matrix(expand.grid(replicate(dim, 0:order, simplify = FALSE)))
}