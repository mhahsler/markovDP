#' Estimate the Convergence Horizon for an Infinite-Horizon MDP
#'
#' Many sampling-based methods require a finite horizon. Fo infinite horizons,
#' discounting leads to convergences during a finite horizon. This function 
#' estimates the number of steps till convergence using rules of thumb.
#'
#' The horizon is estimated differently for the discounted and the undiscounted
#' case.
#' 
#' ## Discounted Case
#'
#' The effect of the largest reward \eqn{R_{\mathrm{max}}} 
#' update decreases with \eqn{t} as \eqn{\delta_t = \gamma^t R_{\mathrm{max}}}.
#' The convergence horizon is estimated as the smallest \eqn{t} for which 
#' \eqn{\delta_t < \delta}. 
#'
#' ## Undiscounted Case
#'
#' For the undiscounted case, episodes end when an absorbing state is reached. 
#' It cannot be guaranteed that a model will reach an absorbing state. 
#' To avoid infinite loops, we set the maximum horizon such that each entry in 
#' the Q-table is on average updated `n_updates` times. This is a very rough 
#' rule ot thumb.
#'
#' @family MDP
#'
#' @param model an MDP model.
#' @param delta maximum update error.
#' @param n_updates integer; average number of time each state is updated. 
#' @return An estimated convergence horizon.
#' @author Michael Hahsler
#' @examples
#' data(Maze)
#' Maze
#'
#' convergence_horizon(Maze)
#'
#' # make the Maze into a discounted problem where future rewards count less.
#' Maze_discounted <- Maze
#' Maze_discounted$discount <- .9
#' Maze_discounted
#'
#' convergence_horizon(Maze_discounted)
#' @export
convergence_horizon <- function(model, delta = 0.001, n_updates = 10) {
  horizon <- model$horizon %||% Inf
  if (!is.infinite(horizon))
    warning("model horizon is not infinite!")
  
  discount <- model$discount %||% 1
  max_abs_R <- max(abs(.reward_range(model)))
  
  if (discount < 1) {
    # steps to get the update to be smaller than delta:
    # delta = gamma^horizon max_abs_R
    horizon <-
      ceiling(log(delta / max_abs_R) / log(model$discount))
  } else {
    # steps to make the influence of max_abs_R for averaging < delta
    horizon <- ceiling(length(S(model)) * length(A(model)) * n_updates)
    warning("discount needs to be <1 to guarantee convergence.\n",
            "  Using a maximum horizon of |S| x |A| x n_updates = ",
            horizon)
  }
 
  horizon 
}
