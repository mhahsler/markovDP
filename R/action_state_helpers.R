#' Conversions for Action and State IDs and Labels
#'
#' Several helper functions to convert state and action (integer) IDs
#' to labels and vice versa.
#' 
#' `normalize_state()` and `normalize_action()` convert labels or ids into a
#' factor which is the standard representation. If only the label or the 
#' integer id (i.e., the index) is needed, the additional functions can be used.
#' These are typically a lot faster.
#' 
#' To support a factored state representation as feature vectors, 
#' `state2features()` and `feature2states()` are provided. 
#' A factored state is represented as a **row** vector for a single 
#' state (ceonveniently created via `s()`) or a matrix with row vectors for a 
#' set of states are used. State labels
#' are constructed in the form `s(feature1, feature2, ...)`. 
#' Factored state representation
#' is used for value function approximation (see [`solve_MDP_APPROX()`]) and
#' for [MDPTF] to describe MDP's via a transition function between factored 
#' states.
#'
#' @name action_state_helpers
#' @aliases action_state_helpers
#'
#' @param state a state labels
#' @param action a action labels
#' @param model a MDP model
#'
#' @returns
#' Functions ending in
#'
#' * `_id` return an integer id,
#' * `_label` return a character string,
#' * `_features` return a state feature matrix,
#' * no ending return a factor.
#'
#' Other functions: 
#' * `state2features()` returns a feature vector/matrix.
#' * `features2state(x)` returns a state label in the format `s(feature list)`.
#' * `s()` returns a state features row vector.
NULL

#' @rdname action_state_helpers
#' @export
normalize_state <- function(state, model) {
  if (is.null(state) || is.factor(state))
    return(state)
  
  # from features
  if (is.matrix(state))
    state <- features2state(state)
  
  if (is.numeric(state)) {
    s <- factor(state, levels = seq_along(S(model)), labels = S(model))
    if (any(is.na(s)))
      stop("Unknown state ", paste(sQuote(state[is.na(s)]), collapse = ", "))
    return(s)
  }
  
  if (is.character(state)) {
    s <- factor(state, levels = S(model))
    if (any(is.na(s)))
      stop("Unknown state ", paste(sQuote(state[is.na(s)]), collapse = ", "))
    return(s)
  }
  stop("Unknown state label: ", sQuote(state))
}

#' @rdname action_state_helpers
#' @export
normalize_state_id <- function(state, model) {
  if (is.null(state) || is.integer(state))
    return(state)
  
  # from features
  if (is.matrix(state))
    state <- features2state(state)
  
  s <- state
  if (is.character(s))
    s <- match(s, S(model))
  
  if (!is.integer(s))
    s <- as.integer(s)
  
  problem <- is.na(s) | s < 1L | s > length(S(model))
  if (any(problem))
    stop("Unknown state ", paste(sQuote(state[problem]), collapse = ", "))
  
  return(s)
}

#' @rdname action_state_helpers
#' @export
normalize_state_label <- function(state, model) {
  if (is.null(state) || is.character(state))
    return(state)
  
  # from features
  if (is.matrix(state))
    return(features2state(state))
  
  if (is.factor(state))
    return(as.character(state))
  
  return(S(model)[state])
}

#' @rdname action_state_helpers
#' @export
normalize_action <- function(action, model) {
  if (is.null(action) || is.factor(action))
    return(action)
  
  if (is.logical(action)) {
    if (length(action) != length(A(model)))
      stop("Illegal action definition (logical)")
    action <- which(action)
  }
  
  if (is.numeric(action)) {
    a <- factor(action, levels = seq_along(A(model)), labels = A(model))
    if (any(is.na(a)))
      stop("Unknown action ", paste(sQuote(action[is.na(a)]), collapse = ", "))
    return(a)
  }
  
  if (is.character(action)) {
    a <- factor(action, levels = A(model))
    if (any(is.na(a)))
      stop("Unknown action ", paste(sQuote(action[is.na(a)]), collapse = ", "))
    return(a)
  }
  
  if (is.null(action))
    return(NULL)
  
  stop("Unknown action label: ", sQuote(action))
}

#' @rdname action_state_helpers
#' @export
normalize_action_id <- function(action, model) {
  if (is.null(action) || is.integer(action))
    return(action)
  
  a <- action
  if (is.character(a))
    a <- match(a, A(model))
  
  if (!is.integer(a))
    a <- as.integer(a)
  
  problem <- is.na(a) | a < 1L | a > length(A(model))
  if (any(problem))
    stop("Unknown action ", paste(sQuote(action[problem]), collapse = ", "))
  
  return(a)
}

#' @rdname action_state_helpers
#' @export
normalize_action_label <- function(action, model) {
  if (is.null(action) || is.character(action))
    return(action)
  
  if (is.factor(action))
    return(as.character(action))
  
  return(A(model)[action])
}

#' @rdname action_state_helpers
#' @export
state2features <- function(state) {
  if (!is.character(state))
    stop("state has to be a state label!")
  
  spl <- strsplit(state, "\\(|,|\\)")
  x <- t(sapply(spl, FUN = function(s) {
    s <- as.numeric(s[-1])
  }))
  
  if (!is.matrix(x) || any(is.na(x)))
    stop("state label ", sQuote(state), 
         " is not formated as 's(feature1, feature2, ...). Cannot extract features!'")
  
  rownames(x) <- state
  
  x
}

#' @rdname action_state_helpers
#' @param x a state feature vector or a matrix of state feature vectors as rows.
#' @export
features2state <- function(x) {
  if (!is.matrix(x))
    stop("Factored state representation has to be a matrix or a row vector!")
      
  paste0("s(", x[, 1], ",", x[, 2], ")")
}

#' @rdname action_state_helpers
#' @param ... features that should be converted into a row
#'            vector used to describe a state.
#' @export
s <- function(...) rbind(c(...))


#' @rdname action_state_helpers
#' @export
normalize_state_features <- function(state, model = NULL) {
  if (is.matrix(state))
    return(state)
  
  if (!is.null(model$state_features)) 
    return(model$state_features[state, ,drop = FALSE])
  
  if (is.factor(state))
    state <- as.character(state)
  
  if (is.character(state))
    return(state2features(state))

  # integer vector is a state ids!
  if (!is.null(S(model)))
    return(state2features(S(model)[state]))
    
  stop("Illegal state specification.")
  }