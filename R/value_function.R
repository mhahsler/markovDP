#' Value Function
#'
#' Extracts the value function from a solved MDP.
#'
#' @family policy
#' @family MDP
#'
#' @param model a solved [MDP].
#' @param drop logical; drop the list for converged, epoch-independent value functions.
#' @param ... further arguments are passed on to [graphics::barplot()]`.
#'
#' @returns the function as a numeric vector with one value for each state.
#'
#' @author Michael Hahsler
#' @keywords hplot
#' @examples
#' data("Maze")
#' sol <- solve_MDP(Maze)
#' sol
#'
#' value_function(sol)
#' plot_value_function(sol)
#'
#' ## finite-horizon problem
#' sol <- solve_MDP(Maze, horizon = 3)
#' policy(sol)
#' value_function(sol)
#' plot_value_function(sol, epoch = 1)
#' plot_value_function(sol, epoch = 2)
#' plot_value_function(sol, epoch = 3)
#'
#' # For a gridworld we can also plot is like this
#' gw_plot(sol, epoch = 1)
#' gw_plot(sol, epoch = 2)
#' gw_plot(sol, epoch = 3)
#' @importFrom graphics barplot text legend
#' @export
value_function <- function(model, drop = TRUE) {
  UseMethod("value_function")
}

#' @export
value_function.MDP <- function(model, drop = TRUE) {
  is_solved_MDP(model, stop = TRUE)
  val <- lapply(policy(model, drop = FALSE), "[[", "U")

  if (drop && length(val) == 1L) {
    val <- val[[1]]
  }

  val
}

#' @rdname value_function
#' @export
#' @param epoch epoch for finite time horizon solutions.
#' @param legend logical; show legend.
#' @param col,ylab,las are passed on to [graphics::barplot()].
#' @param main a main title for the plot. Defaults to the name of the problem.
plot_value_function <- function(model,
                                epoch = 1,
                                legend = TRUE,
                                col = NULL,
                                ylab = "Value",
                                las = 3,
                                main = NULL,
                                ...) {
  UseMethod("plot_value_function")
}

#' @export
plot_value_function.MDP <-
  function(model,
           epoch = 1,
           legend = TRUE,
           col = NULL,
           ylab = "Value",
           las = 3,
           main = NULL,
           ...) {
    is_solved_MDP(model, stop = TRUE)

    if (is.null(main)) {
      main <- paste("Value Function:", model$name, paste0("(", model$solution$method, ")"))
    }

    policy <- policy(model, drop = FALSE)[[epoch]]

    mid <- barplot(
      policy$U,
      col = col,
      ylab = ylab,
      xlab = "State",
      names.arg = paste(model$states),
      las = las,
      main = main,
      ...
    )

    if (legend) {
      text(
        x = mid,
        y = 0,
        labels = policy$action,
        srt = 90,
        adj = c(-.1, .5)
      )
    }
    invisible(NULL)
  }
