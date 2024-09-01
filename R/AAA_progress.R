my_progress_bar <- function(N) 
  progress::progress_bar$new(
    format = "[:bar] :percent (remaining :eta)", 
    total = N)

my_progress_spinner <- function() 
  progress::progress_bar$new(
    format = "(:spin) ticks: :current | elapsed: :elapsed :details", 
    total = NA)