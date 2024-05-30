#' Define an MDP Problem
#'
#' Defines all the elements of a discrete-time finite state-space MDP problem.
#'
#' Markov decision processes (MDPs) are discrete-time stochastic control
#' process. We implement here MDPs with a finite state space. 
#' `MDP()` defines all the element of a MDP problem including the discount rate, the
#' set of states, the set of actions,the transition
#' probabilities, the observation probabilities, and the rewards.
#'
#' In the following we use the following notation. The MDP is a 5-duple:
#'
#' \eqn{(S,A,T,R, \gamma)}.
#'
#' \eqn{S} is the set of states; \eqn{A}
#' is the set of actions; \eqn{T} are the conditional transition probabilities
#' between states; \eqn{R} is the reward function; \eqn{\Omega} is the set of
#' observations; and
#' \eqn{\gamma} is the discount factor. We will use lower case letters to
#' represent a member of a set, e.g., \eqn{s} is a specific state. To refer to
#' the size of a set we will use cardinality, e.g., the number of actions is
#' \eqn{|A|}.
#'
#' **Names used for mathematical symbols in code**
#'
#' * \eqn{S, s, s'}: `'states', start.state', 'end.state'`
#' * \eqn{A, a}: `'actions', 'action'`
#'
#' State names and actions can be specified as strings or index numbers
#' (e.g., `start.state` can be specified as the index of the state in `states`).
#' For the specification as data.frames below, `NA` can be used to mean
#' any  `start.state`, `end.state` or `action`.
#' 
#' **Specification of transition probabilities: \eqn{T(s' | s, a)}**
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
#' * A function with the same arguments are `T_()`, but no default values
#'   that returns the transition probability.
#'
#' **Specification of the reward function: \eqn{R(a, s, s')}**
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
#' * A list of lists. The list levels are `'action'` and `'start.state'`. The list elements
#'   are matrices with
#'   rows representing end states \eqn{s'} and columns representing an observation \eqn{o}.
#'
#' * A function with the same arguments are `R_()`, but no default values
#'   that returns the reward.
#'   
#' To avoid overflow problems with rewards, reward values should stay well within the 
#' range of 
#' `[-1e10, +1e10]`. `-Inf` can be used as the reward for unavailable actions and
#' will be translated into a large negative reward for solvers that only support
#' finite reward values. 
#'
#' Note: The code also includes in `R_()` an argument called `observation`. 
#' Observations are only used POMDPs implemented in package
#' `pomdp` abs must always be `NA` for MDPs.
#'
#' **Start State**
#'
#' The start state of the agent can be a single state or a distribution over the states. 
#' The start state definition is used as the default when the reward is calculated by [reward()]
#' and for simulations with [simulate_MDP()].
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
#' **Time-dependent MDPs**
#'
#' Time dependence of transition probabilities and
#' reward structure can be modeled by considering a set of **episodes**
#' representing **epoch** with the same settings. The length of each episode is
#' specified as a vector for `horizon`, where the length is the number of
#' episodes and each value is the length of the episode in epochs. Transition
#' probabilities, observation probabilities and/or reward structure can contain
#' a list with the values for each episode. The helper function `epoch_to_episode()` converts
#' an epoch to the episode it belongs to.
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
#' @param x a `MDP` object.
#'
#' @return The function returns an object of class MDP which is list with
#'   the model specification. [solve_MDP()] reads the object and adds a list element called
#' `'solution'`.
#' @author Michael Hahsler
#' @examples
#' # Michael's Sleepy Tiger Problem is like the POMDP Tiger problem, but
#' # has completely observable states because the tiger is sleeping in front
#' # of the door. This makes the problem an MDP.
#'
#' STiger <- MDP(
#'   name = "Michael's Sleepy Tiger Problem",
#'   discount = .9,
#'   states = c("tiger-left", "tiger-right"),
#'   actions = c("open-left", "open-right", "do-nothing"),
#'   start = "uniform",
#'
#'   # opening a door resets the problem
#'   transition_prob = list(
#'     "open-left" =  "uniform",
#'     "open-right" = "uniform",
#'     "do-nothing" = "identity"
#'   ),
#'
#'   # the reward helper R_() expects: action, start.state, end.state, value
#'   reward = rbind(
#'     R_("open-left", "tiger-left", v = -100),
#'     R_("open-left", "tiger-right", v = 10),
#'     R_("open-right", "tiger-left", v = 10),
#'     R_("open-right", "tiger-right", v = -100),
#'     R_("do-nothing", v = 0)
#'   )
#' )
#'
#' STiger
#'
#' sol <- solve_MDP(STiger)
#' sol
#'
#' policy(sol)
#' plot_value_function(sol)
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
  # MDP does not have observations
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
    stop("x needs to be a MDP object!")
  }
  solved <- !is.null(x$solution)
  if (stop && !solved) {
    stop("x needs to be a solved MDP. Use solve_MDP() first.")
  }
  
  solved
}

## this is .get_pg_index for MDPs
.get_pol_index <- function(model, epoch) {
  epoch <- as.integer(epoch)
  if (epoch < 1L) {
    stop("Epoch has to be >= 1")
  }
  
  ### (converged) infinite horizon POMDPs. We ignore epoch.
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
  if (any(is.na(id)))
    stop("Unknown state(s): ", paste(state[is.na(id)], collapse = ", "))

  id  
}
  
.get_action_id <- function(model, action) {
  id <- structure(seq_along(model$action), names = model$action)[action]
  if (any(is.na(id)))
    stop("Unknown action(s): ", paste(action[is.na(id)], collapse = ", "))

  id
}

#' @rdname MDP
#' @export
is_timedependent <- function(x) {
  !is.null(x$horizon) && length(x$horizon) > 1L
}


# is a field time-dependent? For time-dependence we have a list of
# matrices/data.frames or for observation_prob we have a list of a list
.is_timedependent_field <- function(x, field) {
  field <-
    match.arg(field, c("transition_prob", "observation_prob", "reward"))
  m <- x[[field]]
  if (is.null(m)) {
    stop("Field ", field, " does not exist.")
  }
  
  # it is a list. time dependent is a list (episodes) of lists
  if (!is.list(m) || is.data.frame(m)) {
    return(FALSE)
  }
  if (!is.list(m[[1]])) {
    return(FALSE)
  }
  
  # time dependent reward is a list of lists of lists
  if (field == "reward" && !is.list(m[[1]][[1]])) {
    return(FALSE)
  }
  
  if (length(m) != length(x$horizon)) {
    stop(
      "Inconsistent POMDP specification. Field ",
      field,
      " does not contain data for the appropriate number of episodes."
    )
  }
  
  TRUE
}

#' @rdname MDP
#' @param epoch integer; an epoch that should be converted to the
#'              corresponding episode in a time-dependent MDP.
#' @export
epoch_to_episode <- function(x, epoch) {
  if (is.null(epoch)) {
    return(1L)
  }
  
  episode <- which(epoch <= cumsum(x$horizon))[1]
  if (is.na(episode)) {
    stop("Epoch does not exist")
  }
  
  episode
}


#' @rdname MDP
#' @param action action as a action label or integer. The value `NA` matches any action.
#' @param start.state,end.state state as a state label or an integer. The value `NA` matches any state.
#' @param observation unused for MDPs. Must be `NA`.
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
           observation = NA,
           value) {
    data.frame(
      action = action,
      start.state = start.state,
      end.state = end.state,
      observation = observation,
      value = as.numeric(value),
      stringsAsFactors = FALSE
    )
  }
