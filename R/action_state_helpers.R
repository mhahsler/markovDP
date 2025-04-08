#' Conversions for Action and State IDs and Labels
#'
#' Several helper functions to convert state and action (integer) IDs
#' to labels and vice versa.
#' 
#' `normalize_state()` and `normalize_action()` convert labels or ids into a
#' desired standard representation. If only the label or the 
#' integer id (i.e., the index) is needed, the additional functions can be used.
#' These are typically a lot faster.
#' 
#' To support a factored state representation as feature vectors, 
#' `state2features()`, `feature2states()`, and `get_state_features()` are provided.
#'  `get_state_features()` is only available if the model explicitly 
#'  stores a finite state space.
#' 
#' **Note:** A factored state is represented as a **row** vector (matrix with a 
#' single row) for a single 
#' state (conveniently created via `s()`) or a matrix with row vectors for a 
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
#' @param state a state in any format
#' @param action a action in any format
#' @param model a MDP model
#' @param as character; specifies the desired output format
#'
#' @returns
#' Functions ending in
#'
#' * `_factor` return a factor,
#' * `_id` return an integer id,
#' * `_label` return a character string,
#' * `_features` return a state feature matrix,
#' * no ending return the type specified as parameter `as`.
#'
#' Other functions: 
#' * `state2features()` returns a feature vector/matrix.
#' * `features2state(x)` returns a state label in the format `s(feature list)`.
#' * `s()` returns a state features row vector.
#' 
#' @examples
#' data(Maze)
#' 
#' # states
#' normalize_state(1, Maze)
#' normalize_state(1, Maze, as = "id")
#' normalize_state(1, Maze, as = "label")
#' normalize_state(1, Maze, as = "features")
#'
#' get_state_features(Maze)
#'
#' # actions
#' normalize_action(1, Maze)
#' normalize_action(1, Maze, as = "id")
#' normalize_action(1, Maze, as = "label")
#'
#' # state label to feature conversion
#' state2features("s(1,1)")
#' s(1,1)
NULL

#' @rdname action_state_helpers
#' @export
normalize_state <- function(state, model, as = "factor") {
  as <- match.arg(as, c("factor", "id", "label", "features"))
  
  switch(as, 
         factor = normalize_state_factor(state, model),
         id = normalize_state_id(state, model),
         label = normalize_state_label(state, model),
         features = normalize_state_features(state, model)
         )
}
  
normalize_state_factor <- function(state, model) {  
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
normalize_state_features <- function(state, model = NULL) {
  if (is.matrix(state)) {
    if (is.null(colnames(state)))
      colnames(state) <- paste0("x", 1:ncol(state))
    return(state)
  }
  
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


#' @rdname action_state_helpers
#' @export
normalize_action <- function(action, model, as = "factor") {
  as <- match.arg(as, c("factor", "id", "label"))
  
  switch(as, 
         factor = normalize_action_factor(action, model),
         id = normalize_action_id(action, model),
         label = normalize_action_label(action, model)
  )
}

normalize_action_factor <- function(action, model) {
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
  colnames(x) <- paste0("x", seq_len(ncol(x)))
  
  x
}

#' @rdname action_state_helpers
#' @param x a state feature vector or a matrix of state feature vectors as rows.
#' @export
features2state <- function(x) {
  if (!is.matrix(x))
    stop("Factored state representation has to be a matrix or a row vector!")
      
  paste0("s(", apply(x, MARGIN = 1, paste0, collapse = ","), ")")
}

#' @rdname action_state_helpers
#' @param ... features that should be converted into a row
#'            vector used to describe a state.
#' @export
s <- function(...) rbind(c(...))

#' @rdname action_state_helpers
#' @export
get_state_features <- function(model) {
  if (!is.null(model$state_features))
    return(model$state_features)
  
  state_features <- NULL
  try(state_features <- state2features(S(model)), silent = TRUE)
  
  return(state_features)
}


# internal: used by solve_MDP_APPROX's transformation functions
get_state_feature_range <- function(model, min = NULL, max = NULL) {
  state_features <- get_state_features(model)
    
  if (!is.null(state_features)) {
    rng <- apply(state_features, MARGIN = 2, range)
    rownames(rng) <- c("min", "max")
    return(rng)
  } else if (model$info$gridworld) {
    # MDPTF can have no state space, check if is a gridworld
    return(rbind(min = 1, max = model$info$dim))
  }
  
  # we have a MDPTF that is not a gridworld at this point
  if (is.null(min) || is.null(max))
    stop("min and max needs to be specified to scale state features to [0,1].")
  
  example_state <- start(model, as = "features")[1, , drop = FALSE]
  if (length(min) == 1)
    min <- rep(min, times = ncol(example_state))
  if (length(max) == 1)
    max <- rep(max, times = ncol(example_state))
  
  if(length(max) != ncol(example_state) ||
     length(min) != ncol(example_state))
    stop("min and max need to be the same length as the number of state features!")
  
  return(rbind(min, max))
}

