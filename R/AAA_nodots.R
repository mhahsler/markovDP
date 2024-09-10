## helper to parse parameter lists with defaults
.nodots <- function(...) {
  l <- list(...)
  if (length(l) > 0L)
    warning("Unknown arguments: ",
            paste(names(l), "=", l, collapse = ", "))
}