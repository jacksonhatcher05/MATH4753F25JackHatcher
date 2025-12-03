#' Central Limit Theorem for Binomial Distribution
#'
#' This function simulates the Central Limit Theorem for a Binomial distribution.
#' It generates sample means and plots a histogram with a theoretical normal curve.
#'
#' @param n Integer. Sample size for each iteration.
#' @param iter Integer. Number of iterations/samples to simulate.
#' @param p Numeric. Probability of success in the Binomial distribution.
#' @param ... Additional arguments passed to hist().
#'
#' @return None. Produces a plot of the sample means and theoretical normal curve.
#' @examples
#' mycltb(n=5, iter=1000, p=0.3)
#' @export
mycltb <- function(n, iter, p=0.5, ...) {
  y <- rbinom(n*iter, size=n, prob=p)
  data <- matrix(y, nr=n, nc=iter, byrow=TRUE)
  w <- apply(data, 2, mean)
  param <- hist(w, plot=FALSE)
  ymax <- 1.1 * max(param$density)

  hist(w, freq=FALSE, ylim=c(0, ymax),
       main=paste("Histogram of sample mean\nn =", n, ", p =", p),
       xlab="Sample mean", ...)

  curve(dnorm(x, mean=n*p, sd=sqrt(n*p*(1-p))), add=TRUE,
        col="Red", lty=2, lwd=3)
}

