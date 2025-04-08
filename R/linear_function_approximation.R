# Linear function approximation for v, q

# all approximations are lists with:
# * x(s, a) .... function to construct from state features, action-state features.
# * f(s, a, w) ... approx. function
# * gradient(s, a, w) ... gradient of f at w
# * w_init ... initial weight vector (typically all 0s)
# * transformation ... a transformation kernel function that is applied to state features in x.
#
#  Note for v: a is null and ignored for approximating v

#' @rdname solve_MDP_APPROX
#' @param state_features a matrix with state features. Each row is the feature
#'    vector for a state.
#' @param transformation a transformation function that is applied to the feature vector
#'    before it is used in the linear approximation. See [transformation].
#' @export
add_linear_approx_Q_function <- function(model,
                                         state_features = NULL,
                                         transformation = transformation_linear_basis,
                                         ...) {
  
  if (!is.null(state_features)) {
    if (!is.matrix(state_features) ||
        !nrow(state_features) == c(length(S(model))))
      stop("state_features needs to be a matrix with one row per state!")
    
    model$state_featues <- state_features
  }
  
  model$approx_Q_function <- linear_function_approx(model, "q", transformation, ...)
  model
}

# internal
linear_function_approx <- function(model,
                                   what = "q",
                                   transformation = transformation_linear_basis,
                                   ...) {
  what <- match.arg(what, c("q", "v"))
  
  transformation <- transformation(model, ...)
  
  state_template <- transformation(start(model, as = "features")[1, ])
  dim_s <- length(state_template)
  
  if (what == "q") {
    n_A <- length(A(model))
    
    # s are state features. Convert to action-state feature.
    x <- function(s, a) {
      # use a transformation function
      s <- transformation(s)
      a <- normalize_action_id(a, model)
      
      # one component per action
      x <- numeric(n_A * dim_s)
      
      a_pos <- 1L + (a - 1L) * dim_s
      x[a_pos:(a_pos + dim_s - 1L)] <- s
      
      x
    }
    
    w_init <- setNames(numeric(n_A * dim_s), paste(
      rep(A(model), each = dim_s), names(state_template), sep = "."))
    
  } else if (what == "v") {
    # a is ignored!
    x <- function(s, a) {
      transformation(s)
    }
    
    w_init <- setNames(numeric(dim_s), names(state_template))
    
  } else stop("Unknown approximation target.")
  
  list(
    transformation = transformation,
    w_init = w_init,
    f = function(s, a, w)
      sum(w * x(s, a)),
    gradient = function(s, a, w)
      x(s, a),
     x = x
  )
  
}

# internal: f can be directly called in code
approx_value <- function(f,
                         state,
                         action = NULL,   # NULL is for v approx.
                         w) {
  
   f$f(state, action, w)
}



#' @rdname solve_MDP_APPROX
#' @param state a state (index or name)
#' @param action an action (index or name)
#' @param w a weight vector
#' @export
approx_Q_value <- function(model,
                           state = NULL,
                           action = NULL,
                           w = NULL) {
  action <- action %||% A(model)
  
  state <- state %||% S(model)
  if (is.null(state))
    stop("a state needs to be specified!")
  
  w <- w %||% model$solution$w
  if (is.null(w))
    stop("weight vector w is missing.")
  
  if(!is.matrix(state))
    state <- normalize_state_features(state, model)
  
  if (length(action) > 1L)
    return(sapply(
      action,
      FUN = function(a)
        approx_Q_value(model, state, a, w)
    ))
  
  if (nrow(state) == 1L)
    model$approx_Q_function$f(drop(state), action, w)
  else
    apply(state, MARGIN = 1, model$approx_Q_function$f, action, w)
}

#' @rdname solve_MDP_APPROX
#' @export
approx_greedy_action <- function(model,
                                 state,
                                 w = NULL,
                                 epsilon = 0) {
  if (epsilon == 0 || runif(1) > epsilon) {
    a <- which.max.random(approx_Q_value(model, state, w = w))
  } else {
    a <- sample.int(length(A(model)), 1L)
  }
  
  normalize_action(a, model)
}

#' @rdname solve_MDP_APPROX
#' @export
approx_greedy_policy <- function(model, w = NULL) {
  w <- w %||% model$solution$w
  S <- S(model)
  
  if (is.null(S))
    stop("Policy is only available for a specified finite state space!")
  
  qs <- approx_Q_value(model, w = w)
  
  data.frame(
    state = S,
    V = apply(qs, MARGIN = 1, max),
    action = normalize_action(apply(qs, MARGIN = 1, which.max.random), model),
    row.names = NULL
  )
}


# internal: transpose and reorder rows for a proper image
#' @importFrom graphics axis contour image
pimage <- function (x1,
                    x2,
                    z,
                    col = hcl.colors(12, "YlOrRd", rev = TRUE),
                    image = TRUE,
                    contour = TRUE,
                    axes = TRUE,
                    ...) {
  
  z <- t(z)[, rev(seq_len(nrow(z))), drop = FALSE]
  
  if (image) {
    image(x2, x1, z, axes = FALSE, col = col, ...)
    
    if (contour)
      contour(x2, x1, z, axes = FALSE, add = TRUE)
  } else
    contour(x2, x1, z, axes = FALSE, ...)
  
  box()
  
  if (axes) {
    p <- pretty(x2)
    axis(1L,
         at = seq(min(p), max(p), along.with = p),
         labels = p)
    p <- pretty(x1)
    axis(2L,
         at = seq(min(p), max(p), along.with = p),
         labels = p)
  }
}

#' @rdname solve_MDP_APPROX
#' @param image,contour logical; include the false color image or the 
#'        contours in the plot?
#' @param main title for the plot.
#' @param res resolution as the number of values sampled from each feature.
#' @param col colors for the passed on to [`image()`].
#' @export
approx_V_plot <- function(model,
                          min = NULL,
                          max = NULL,
                          w = NULL,
                          res = 25,
                          col = hcl.colors(res, "YlOrRd", rev = TRUE),
                          image = TRUE,
                          contour = TRUE,
                          main = NULL,
                          ...) {
  if (is.null(main)) {
    main <- paste("Approx. Value Function:",
                  model$name,
                  paste0("(", model$solution$method, ")"))
  }
  
  rng <- get_state_feature_range(model, min, max)
  min <- rng[1, ]
  max <- rng[2, ]
  dim_s <- length(max)
  
  if (dim_s != 2)
    stop("This visualization only works for states with 2 features!")
  
  x1 <- seq(min[1], max[1], length.out = res)
  x2 <- seq(min[2], max[2], length.out = res)
  
  V_val <- Vectorize(function(x1, x2)
    max(approx_Q_value(model, state = s(x1, x2), w = w)))
  V <- outer(x1, x2, V_val)
  
  pimage(x1, x2, V, col, image, contour, main = main, ...)
}
