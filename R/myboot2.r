#' Bootstrap Confidence Intervals
#'
#' Computes bootstrap confidence intervals for a given statistic and produces a histogram
#' of the bootstrap replicates. The function allows resampling from the input data and
#' calculates point estimates and percentile-based confidence intervals.
#'
#' @param iter Integer. Number of bootstrap iterations (default = 10000).
#' @param x Numeric vector. The input data sample.
#' @param fun Function. Statistic function to compute for each bootstrap sample (default = mean).
#' @param alpha Numeric. Significance level for confidence interval (default = 0.05 for 95% CI).
#' @param cx Numeric. Scaling factor for text size in plots (default = 1.5).
#' @param ... Additional graphical parameters passed to the \code{hist()} function.
#'
#' @return A list containing:
#' \describe{
#'   \item{ci}{Numeric vector. The bootstrap confidence interval for the statistic.}
#'   \item{pte}{Numeric. Point estimate of the statistic for the original sample.}
#'   \item{xstat}{Numeric vector. Bootstrap replicates of the statistic.}
#'   \item{x}{Original input data vector.}
#' }
#'
#' @export
#'
#' @examples
#' set.seed(10)
#' sample_data <- rnorm(20, mean = 5, sd = 2)
#' res <- myboot2(x = sample_data, fun = mean, alpha = 0.05)
#' res$ci   # Bootstrap 95% CI
#' res$pte  # Point estimate
#' hist(res$xstat)  # Histogram of bootstrap replicates
myboot2 <- function(iter = 10000, x, fun = mean, alpha = 0.05, cx = 1.5, ...) {
  n <- length(x)
  y <- sample(x, n * iter, replace = TRUE)

  # Use full argument names
  rs.mat <- matrix(y, nrow = n, ncol = iter, byrow = TRUE)
  xstat <- apply(rs.mat, 2, fun)

  ci <- quantile(xstat, c(alpha/2, 1 - alpha/2))

  para <- hist(xstat, freq = FALSE, las = 1,
               main = paste("Histogram of Bootstrap sample statistics",
                            "\n", "alpha=", alpha, " iter=", iter, sep = ""),
               ...)

  mat <- matrix(x, nrow = length(x), ncol = 1, byrow = TRUE)
  pte <- apply(mat, 2, fun)
  abline(v = pte, lwd = 3, col = "Black")
  segments(ci[1], 0, ci[2], 0, lwd = 4)
  text(ci[1], 0, paste("(", round(ci[1], 2), sep = ""), col = "Red", cex = cx)
  text(ci[2], 0, paste(round(ci[2], 2), ")", sep = ""), col = "Red", cex = cx)
  text(pte, max(para$density) / 2, round(pte, 2), cex = cx)

  invisible(list(ci = ci, pte = pte, xstat = xstat, x = x))
}
