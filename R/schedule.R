#' Schedules to Reduce Alpha, Epsilon and Other Parameters
#'
#' Several schedule functions to reduce learning parameters are available using 
#' generator functions.
#'
#' @family solver
#' @name schedule
#' @aliases schedule
#' 
#' @details
#' Several learning parameters need to be reduced during learning to ensure
#' convergence. We provide several schedule function generators 
#' that reduce learning 
#' parameters after each processed episode or according to the number of 
#' times a state-action combination was tried.
#' 
#' Here are the definitions of the available schedules:
#' 
#' ```{r comment = "", results = "asis", echo = FALSE} 
#' schedules <- list(
#'    schedule_exp = schedule_exp, 
#'    schedule_exp2 = schedule_exp2,
#'    schedule_log = schedule_log,
#'    schedule_linear = schedule_linear,
#'    schedule_harmonic = schedule_harmonic
#'    ) 
#'
#' for(i in seq_along(schedules)) {
#'  cat("\n", "* ", names(schedules)[i], ": ")
#'  cat("`", as.character(body(schedules[[i]])[-1]), "`")
#'  cat("\n")
#' }
#' ```
#' 
#' `t` is the time step, episode number or count starting with 1.
#' 
#' @param start start value for the schedule.
#' @param end end value for the schedule.
#' @param decay decay factor for exponential schedules.
#' @param basis basis for exponential schedules.
#' @param n number of steps (e.g., epochs) for the schedule.
#' 
#' @examples
#' # create an exponential  schedule function
#' s_exp <- schedule_exp(1, decay = .1)
#' s_exp
#' 
#' # plot the schedule for 100 episodes.
#' episode <- seq_len(100)
#' plot(x = episode, s_exp(1:100), type = "l")
#' 
#' # compare some schedule examples
#' schedules <- cbind(
#'   `exp decay = 0.1` = schedule_exp(1, decay = .1)(1:100),
#'   `exp decay = 0.01` = schedule_exp(1, decay = .01)(1:100),
#'   `exp basis = 0.9` = schedule_exp2(1, basis = .9)(1:100),
#'   `linear` = schedule_linear(1, end = 0, n = 100)(1:100),
#'   `log` = schedule_log(1)(1:100),
#'   `harmonic` = schedule_harmonic(1, n = 100)(1:100),
#'   `harmonic start = 10` = schedule_harmonic(10, n = 100)(1:100)
#'   )
#'  
#' matplot(schedules, type = "l", 
#'         col = 1:ncol(schedules),
#'         lty = 1:ncol(schedules)
#'         )
#' legend("topright", 
#'         legend = colnames(schedules), 
#'         col = 1:ncol(schedules), 
#'         lty = 1:ncol(schedules),
#'         cex = 0.8
#'         )
#' @export
schedule_exp <- function(start, decay) { 
  function(t) start * exp(-decay * (t-1))
}

#' @rdname schedule
#' @export
schedule_exp2 <- function(start, basis) { 
  function(t) start * basis ^ (t-1)
}

#' @rdname schedule
#' @export
schedule_log <- function(start) { 
  function(t) pmin(start / log(1 + t), start)
}

#' @rdname schedule
#' @export
schedule_linear <- function(start, end, n) {
  step <- (start - end) / n; 
  function(t) start - t * step
}

#' @rdname schedule
#' @export
schedule_harmonic <- function(start, n) {
  function(t) pmin(start / t, 1) 
}


# internal function
show_schedule <- function(schedule) {
  if (!is.function(schedule))
    return (schedule)
  
  schedule_env <- as.list(environment(schedule))
  paste0(paste(deparse1(schedule_exp(), collapse = " "), 
                  "with", 
                  paste0(names(schedule_env), " = ", 
                         unlist(schedule_env), collapse = ", ")))
}
