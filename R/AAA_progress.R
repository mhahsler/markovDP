my_progress_bar <- function(N, name = NULL, format_extra = NULL, 
                            clear = FALSE, ...) 
  progress::progress_bar$new(
    format = paste(name, "[:bar] :percent (remaining :eta | elapsed :elapsed)", format_extra),
    clear = clear, 
    total = N,
    ...)

my_progress_spinner <- function(name = NULL, ticks = "iterations", 
                                format_extra = "", clear = FALSE, ...) 
  progress::progress_bar$new(
    format = paste(name, "(:spin)", paste0(ticks, ":"), " :current | elapsed: :elapsed", format_extra), 
    clear = clear,
    total = NA,
    ...)