#' myncurve
#'
#' @param mu Mean of the normal distribution
#' @param sigma Standard deviation of the normal distribution
#' @param a Upper limit for shading
#'
#' @returns A list with mu, sigma, and probability P(X <= a)
#' @export
#'
#' @examples
#' myncurve(mu = 10, sigma = 5, a = 6)
myncurve <- function(mu, sigma, a) {
  # Draw the normal curve
  curve(dnorm(x, mean = mu, sd = sigma),
        xlim = c(mu - 3*sigma, mu + 3*sigma),
        main = paste("Normal(", mu, ",", sigma, ")"),
        ylab = "Density", xlab = "X")

  # Shade the area under the curve from left tail up to a
  x_vals <- seq(mu - 3*sigma, a, length.out = 1000)
  y_vals <- dnorm(x_vals, mean = mu, sd = sigma)
  polygon(c(x_vals, a), c(y_vals, 0), col = rgb(0, 0, 1, 0.3), border = NA)

  # Compute probability P(X <= a)
  prob <- pnorm(a, mean = mu, sd = sigma)

  # Return results as a named list
  return(list(mu = mu, sigma = sigma, probability = prob))
}
