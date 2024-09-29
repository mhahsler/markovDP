#' Define an MDP Problem
#'
#' Defines all the elements of a discrete-time finite state-space MDP problem.
#'
#' Markov decision processes (MDPs) are discrete-time stochastic control
#' process. We implement here MDPs with a finite state space.
#' `MDP()` defines all the element of an MDP problem including the discount rate, the
#' set of states, the set of actions,the transition
#' probabilities, and the rewards.
#'
#' In the following we use the following notation. The MDP is a 5-duple:
#'
#' \eqn{(S,A,T,R, \gamma)}.
#'
#' \eqn{S} is the set of states; \eqn{A}
#' is the set of actions; \eqn{T} are the conditional transition probabilities
#' between states; \eqn{R} is the reward function; and
#' \eqn{\gamma} is the discount factor. We will use lower case letters to
#' represent a member of a set, e.g., \eqn{s} is a specific state. To refer to
#' the size of a set we will use cardinality, e.g., the number of actions is
#' \eqn{|A|}.
#'
#' ## Names used for mathematical symbols in code
#'
#' * \eqn{S, s, s'}: `'states', start.state', 'end.state'`
#' * \eqn{A, a}: `'actions', 'action'`
#'
#' State names and actions can be specified as strings or index numbers
#' (e.g., `start.state` can be specified as the index of the state in `states`).
#' For the specification as data.frames below, `NA` can be used to mean
#' any  `start.state`, `end.state` or `action`.
#'
#' ## Specification of transition model: \eqn{T(s' | s, a)}
#'
#' Transition probability to transition to state \eqn{s'} from given state \eqn{s}
#' and action \eqn{a}. The transition probabilities can be
#' specified in the following ways:
#'
#' * A data.frame with columns exactly like the arguments of `T_()`.
#'   You can use `rbind()` with helper function `T_()` to create this data
#'   frame. Probabilities can be specified multiple times and the definition that
#'   appears last in the data.frame will take affect.
#'
#' * A named list of matrices, one for each action. Each matrix is square with
#'   rows representing start states \eqn{s} and columns representing end states \eqn{s'}.
#'   Instead of a matrix, also the strings `'identity'` or `'uniform'` can be specified.
#'
#' * One of the following functions:
#'    - A function with the argument list `model`, `action`, `start.state`, `end.state` 
#'      which returns a single transition probability.
#'    - A function with the argument list `model`, `action`, `start.state` 
#'      which returns a transition probability vector for all end states. This vector can
#'      be dense, a [Matrix::sparseVector] or a named vector only containing 
#'      the non-zero probabilities named by the corresponding end state.
#'    
#'    The arguments `action`, `start.state`, and `end.state` will be called with the 
#'    state name as a string.
#'
#' ## Specification of the reward function: \eqn{R(a, s, s')}
#'
#' The reward function can be specified in the following
#' ways:
#'
#' * A data frame with columns named exactly like the arguments of `R_()`.
#'   You can use `rbind()`
#'   with helper function `R_()` to create this data frame. Rewards can be specified
#'   multiple times and the definition that
#'   appears last in the data.frame will take affect.
#'
#' * A list of state x state matrices. 
#'   The list elements are for `'action'`. The matrix rows are `start.state` and the 
#'   columns are `end.state`. 
#'
#' * A function following the same rules as for transition probabilities.
#'
#' To avoid overflow problems with rewards, reward values should stay well within the
#' range of
#' `[-1e10, +1e10]`. `-Inf` can be used as the reward for unavailable actions and
#' will be translated into a large negative reward for solvers that only support
#' finite reward values.
#'
#' ## Specification of the Start State
#'
#' The start state of the agent can be a single state or a distribution over the states.
#' The start state definition is used as the default when the reward is calculated by [reward()]
#' and for sampling with [sample_MDP()].
#'
#' Options to specify the start state are:
#'
#' * A string specifying the name of a single starting state.
#' * An integer in the range \eqn{1} to \eqn{n} to specify the index of a single starting state.
#' * The string `"uniform"` where the start state is chosen using a uniform distribution over all states.
#' * A probability distribution over the states. That is, a vector
#'   of \eqn{|S|} probabilities, that add up to \eqn{1}.
#'
#' The default state state is a uniform
#' distribution over all states.
#' 
#' ## Normalizing MDP Descriptions
#' 
#' Different components can be specified in various ways. It is often 
#' necessary to convert each component into a specific form (e.g., a 
#' dense matrix) to save time during access.  
#' Convert the Complete MDP Description into a consistent form
#' `normalize_MDP()` converts all components of the MDP description
#'  into a consistent form and
#' returns a new MDP definition where `transition_prob`,
#' `reward`, and `start` are normalized. This includes the internal
#' representation (dense, sparse, as a data.frame) and 
#' also, `states`, and `actions` are ordered as given in the problem
#' definition to make safe access using numerical indices possible. Normalized
#' MDP descriptions can be
#' used in custom code that expects consistently a certain format.
#'
#' @family MDP
#' @family MDP_examples
#'
#' @param states a character vector specifying the names of the states.
#' @param actions a character vector specifying the names of the available
#' actions.
#' @param transition_prob Specifies the transition probabilities between
#' states.
#' @param reward Specifies the rewards dependent on action and states.
#' @param discount numeric; discount rate between 0 and 1.
#' @param horizon numeric; Number of epochs. `Inf` specifies an infinite
#' horizon.
#' @param start Specifies in which state the MDP starts.
#' @param info A list with additional information.
#' @param name a string to identify the MDP problem.
#' @param model,x a `MDP` object.
#'
#' @return The function returns an object of class MDP which is list with
#'   the model specification. [solve_MDP()] reads the object and adds a list element called
#' `'solution'`.
#' @author Michael Hahsler
#' @examples
#' # simple MDP example
#' #
#' # states:    s1 s2 s3 s4
#' # transitions: forward moves -> and backward moves <-
#' # start: s1
#' # reward: s1, s2, s4 = 0 and s3 = 1
#'
#' car <- MDP(
#'   states = c("s1", "s2", "s3", "s4"),
#'   actions = c("forward", "back", "stop"),
#'   transition <- list(
#'     forward = rbind(c(0, 1, 0, 0), 
#'                     c(0, 0, 1, 0), 
#'                     c(0, 0, 0, 1), 
#'                     c(0, 0, 0, 1)),
#'     back =    rbind(c(1, 0, 0, 0), 
#'                     c(1, 0, 0, 0), 
#'                     c(0, 1, 0, 0), 
#'                     c(0, 0, 1, 0)),
#'     stop = "identity"
#'   ),
#'   reward = rbind(
#'     R_(value = 0),
#'     R_(end.state = "s3", value = 1)
#'   ),
#'   discount = 0.9,
#'   start = "s1",
#'   name = "Simple Car MDP"
#' )
#'
#' car
#'
#' # internal representation
#' str(car)
#'
#' # accessing elements
#' transition_matrix(car, sparse = FALSE)
#' transition_matrix(car, sparse = TRUE)
#' reward_matrix(car, sparse = FALSE)
#' reward_matrix(car, sparse = TRUE)
#'
#' sol <- solve_MDP(car)
#' policy(sol)
#' 
#' # normalize MDP: make everything dense (transition_prob, reward and start)
#' car_dense <- normalize_MDP(car, sparse = FALSE)
#' str(car_dense)
#' @export
MDP <- function(states,
                actions,
                transition_prob,
                reward,
                discount = .9,
                horizon = Inf,
                start = "uniform",
                info = NULL,
                name = NA) {
  # MDP does not have observations (only POMDPs in package pomdp have it)
  if (is.data.frame(reward)) {
    if (!is.null(reward$observation) && !all(is.na(reward$observation))) {
      stop("MDPs do not have observations. Remove observation information from rewards!")
    }
    reward$observation <- NULL
  }

  x <- list(
    name = name,
    discount = discount,
    horizon = horizon,
    states = states,
    actions = actions,
    transition_prob = transition_prob,
    reward = reward,
    info = info,
    start = start
  )

  class(x) <- list("MDP", "list")
  x <- check_and_fix_MDP(x)
    
  # this takes a while
  #x$absorbing_states <- absorbing_states(x)
  #x$unreachable_states <- unreachable_states(x)
  
  x
}


#' @export
print.MDP <- function(x, ...) {
  writeLines(paste(
    paste(class(x), collapse = ", "),
    "-",
    x$name
  ))

  if (!is.null(x$discount)) {
    writeLines(sprintf(
      "  Discount factor: %s",
      paste(x$discount, collapse = "+")
    ))
  }

  if (!is.null(x$horizon)) {
    writeLines(sprintf(
      "  Horizon: %s epochs",
      paste(x$horizon, collapse = " + ")
    ))
  }

  writeLines(sprintf(
    "  Size: %d states / %d actions",
    length(x$states),
    length(x$actions)
  ))

  writeLines(paste0("  Start: ", shorten(paste(
    x$start,
    collapse = ", "
  ), n = -10L)))

  if (is_solved_MDP(x)) {
    writeLines(c(
      "  Solved:",
      sprintf(
        "    Method: %s",
        sQuote(x$solution$method)
      ),
      sprintf(
        "    Solution converged: %s",
        x$solution$converged
      )
    ))
  }

  writeLines("")

  writeLines(strwrap(
    paste("List components:", paste(sQuote(names(
      x
    )), collapse = ", "), "\n"),
    indent = 2,
    exdent = 4
  ))
}


#' @rdname MDP
#' @param stop logical; stop with an error.
#' @export
is_solved_MDP <- function(x, stop = FALSE) {
  if (!inherits(x, "MDP")) {
    stop("x needs to be an MDP object!")
  }
  solved <- !is.null(x$solution)
  if (stop && !solved) {
    stop("x needs to contain a policy. Use solve_MDP() or add_policy() first.")
  }

  solved
}

## this is .get_pg_index for MDPs
.get_pol_index <- function(model, epoch) {
  epoch <- as.integer(epoch)
  if (epoch < 1L) {
    stop("Epoch has to be >= 1")
  }

  ### (converged) infinite horizon MDP. We ignore epoch.
  if (length(model$solution$policy) == 1L) {
    return(1L)
  }

  ### regular epoch for finite/infinite horizon case
  if (epoch > length(model$solution$policy)) {
    stop(
      "MDP model has only a policy up to epoch ",
      length(model$solution$policy)
    )
  }

  return(epoch)
}

# convert names or ids to ids!
.get_state_id <- function(model, state) {
  id <- structure(seq_along(model$states), names = model$states)[state]
  if (any(is.na(id))) {
    stop("Unknown state(s): ", paste(state[is.na(id)], collapse = ", "))
  }

  id
}

.get_action_id <- function(model, action) {
  id <- structure(seq_along(model$action), names = model$action)[action]
  if (any(is.na(id))) {
    stop("Unknown action(s): ", paste(action[is.na(id)], collapse = ", "))
  }

  id
}

#' @rdname MDP
#' @param action action as a action label or integer. The value `NA` matches any action.
#' @param start.state,end.state state as a state label or an integer. The value `NA` matches any state.
#' @param probability,value Values used in the helper functions `T_()` and `R_()`.
#'
#' @export
T_ <-
  function(action = NA,
           start.state = NA,
           end.state = NA,
           probability) {
    data.frame(
      action = action,
      start.state = start.state,
      end.state = end.state,
      probability = as.numeric(probability),
      stringsAsFactors = FALSE
    )
  }

#' @rdname MDP
#' @export
R_ <-
  function(action = NA,
           start.state = NA,
           end.state = NA,
           value) {
    data.frame(
      action = action,
      start.state = start.state,
      end.state = end.state,
      value = as.numeric(value),
      stringsAsFactors = FALSE
    )
  }

#' @rdname MDP
#' @param epoch integer; an epoch that should be converted to the
#'              corresponding episode in a time-dependent MDP.
#' @export
epoch_to_episode <- function(x, epoch) {
  UseMethod("epoch_to_episode")
}

#' @export
epoch_to_episode.MDP <- function(x, epoch) {
  if (is.null(epoch)) {
    return(1L)
  }

  episode <- which(epoch <= cumsum(x$horizon))[1]
  if (is.na(episode)) {
    stop("Epoch does not exist")
  }

  # MDP
  if (episode != 1) {
    stop("MDPs do not support time dependence!")
  }

  episode
}



#' @rdname MDP
#' @param sparse logical; use sparse representation. matrices when the density is below 50% and keeps data.frame representation
#'  for the reward field. 
#' @param trans_function logical; convert functions into matrices?
#' @param trans_keyword logical; convert distribution keywords (uniform and identity)
#'  in `transition_prob` matrices?
#' @param keep_reward_df logical; if reward is a data.frame, then keep it.
#' @param precompute_absorbing_unreachable logical; should absorbing and unreachable states be precalculated?
#' @param progress logical; show a progress bar with estimated time for completion.
#' @export
normalize_MDP <- function(model,
                          sparse = TRUE,
                          trans_keyword = TRUE,
                          trans_function = TRUE,
                          keep_reward_df = FALSE,
                          precompute_absorbing_unreachable = TRUE,
                          progress = TRUE
                          ) {
  if (!inherits(model, "MDP")) {
    stop("model is not an MDP object!")
  }
  
  # start state vector + transitions matrix + reward matrix + check and fix
  n_states <- length(model$states)
  n_actions <- length(model$actions)
  N <- n_states + n_actions * n_states * n_states * 3
  if (progress) {
    pb <- my_progress_bar(N, name = "normalize_MDP")
    pb$tick(0)
  }
  
  if (trans_keyword)
    model$start <- start_vector(model, sparse = sparse)
  
  if (progress)
    pb$tick(n_states)
  
  if (is.function(model$transition_prob) && !trans_function) {
    # do nothing
  } else {
    #model$transition_prob <-
    #  transition_matrix(model, sparse = sparse, trans_keyword = trans_keyword)
    
    model$transition_prob <- 
      sapply(
        model$actions,
        FUN = function(a) {
          tm <- transition_matrix(model, a, sparse = sparse)
          if (progress)
            pb$tick(n_states * n_states)
          tm
        },
        simplify = FALSE,
        USE.NAMES = TRUE
      )
  }
  
  #  if (progress)
  #    pb$tick(n_actions * n_states * n_states)
  
  if ((is.function(model$reward) && !trans_function)) {
    # do nothing
  } else if (is.data.frame(model$reward) && keep_reward_df) {
    # do nothing
  } else {
    model$reward <- reward_matrix(model, sparse = sparse)
  }
  
  if (progress)
    pb$tick(n_actions * n_states * n_states)
  
  # make sure order is OK
  model <- check_and_fix_MDP(model)
  
  if (progress)
    pb$tick(n_actions * n_states * n_states)
  
  # TODO: remember absorbing states
  # remember recalculated absorbing/unreachable states
  model$absorbing_states <- NULL
  model$unreachable_states <- NULL
  if (precompute_absorbing_unreachable) {
    model$absorbing_states <- absorbing_states(model, sparse = TRUE)
    model$unreachable_states <- unreachable_states(model, sparse = TRUE)
  }
  model
}
