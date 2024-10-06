my_progress_bar <- function(N,
                            name = NULL,
                            format_extra = NULL,
                            clear = FALSE,
                            show_after = 0.5,
                            ...)
  progress::progress_bar$new(
    format = paste(
      name,
      "[:bar] :percent (remaining :eta | elapsed :elapsed)",
      format_extra
    ),
    total = N,
    clear = clear,
    show_after = show_after,
    ...
  )

my_progress_spinner <- function(name = NULL,
                                ticks = "iterations",
                                format_extra = "",
                                clear = FALSE,
                                show_after = 0.5,
                                ...)
  progress::progress_bar$new(
    format = paste(
      name,
      "(:spin)",
      paste0(ticks, ":"),
      " :current | elapsed: :elapsed",
      format_extra
    ),
    total = NA,
    clear = clear,
    show_after = show_after,
    ...
  )