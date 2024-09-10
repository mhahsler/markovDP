# Accessor Functions for transitions and observations
#
# Representations:
# Default:
# * Sparse (list):
#     Trans: A action list -> start.state x end.state sparse matrix
#
# Others
# * Dense (list): Same as sparse with dense matrices
# * df: A data.frame with value
# * A function can be converted to a list
#
# sparse = NULL translates functions/data frames/strings
#

#' @include accessors.R
#' @rdname accessors
#' @param trans_keyword logical; translate keywords like "uniform" into matrices.
#' @export
transition_matrix <- function(model,
                              action = NULL,
                              start.state = NULL,
                              end.state = NULL,
                              ...,
                              sparse = NULL,
                              trans_keyword = TRUE) {
  UseMethod("transition_matrix")
}

#' @export
transition_matrix.MDP <-
  function(model,
           action = NULL,
           start.state = NULL,
           end.state = NULL,
           ...,
           sparse = NULL,
           trans_keyword = TRUE) {
    value_matrix(
      model,
      "transition_prob",
      action,
      start.state,
      end.state,
      sparse,
      trans_keyword
    )
  }
