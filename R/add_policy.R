#' Add a Policy to a MDP Problem Description
#'
#' Add a policy to a MDP problem description allows the user to
#' test policies on modified problem descriptions or to test manually created
#' policies.
#'
#' The new policy needs to be a data.frame with one row for each state in the
#' order the states are defined in the model. The only required column is
#'
#' * `action`: the action prescribed in the state corresponding to the row.
#'
#' Optional columns are
#'
#' * `state`: the state names in the order of the states in the model.
#'    The needed names can be obtained by from the `$states` element of the model.
#' * `U`: with the utility given by the value function for the state.
#'
#' @family MDP
#'
#' @param model a [MDP] model description.
#' @param policy a policy data.frame.
#'
#' @return The model description with the added policy.
#'
#' @author Michael Hahsler
#' @examples
#' data(Maze)
#'
#' sol <- solve_MDP(Maze)
#' sol
#'
#' policy(sol)
#' reward(sol)
#'
#' # Add a random policy
#' random_pol <- random_policy(Maze)
#' random_pol
#' sol_random <- add_policy(Maze, random_pol)
#' policy(sol_random)
#' reward(sol_random)
#' @export
add_policy <- function(model, policy) {
  UseMethod("add_policy")
}

#' @export
add_policy.MDP <- function(model, policy) {
  if (inherits(policy, "MDP")) {
    policy <- policy(policy)
  }

  if (is.null(policy$U)) {
    #policy$U <- policy_evaluation(model, policy)
    policy$U <- NA
  }

  solution <- list(
    method = "manual",
    policy = list(policy),
    converged = NA
  )

  model$solution <- solution
  model <- check_and_fix_MDP(model)

  model
}
