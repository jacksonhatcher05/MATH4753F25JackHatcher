#' Birthday Function
#'
#' @param x Integer or Vector, Number of people in group.
#'
#' @returns Values between 0 and 1 giving the Probability
#' @export
#'
#' @examples
#' birthday(20)
#' birthday(20:24)
#'
birthday <- function(x) {
  1 - exp(lchoose(365, x) + lfactorial(x) - x * log(365))
}

birthday(20:25)
