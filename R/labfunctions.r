# R/lab5functions.R

#' Simulate a Binomial Distribution
#'
#' This function simulates a binomial experiment for a given number of iterations,
#' sample size, and probability of success.
#'
#' @param iter Number of simulations to run
#' @param n Number of trials per simulation
#' @param p Probability of success for each trial
#'
#' @return A vector of relative frequencies for 0:n successes
#' @export
#'
#' @examples
#' mybin(iter = 100, n = 10, p = 0.7)
mybin <- function(iter = 100, n = 10, p = 0.5) {
  results <- numeric(n + 1)
  for (i in 1:iter) {
    x <- rbinom(1, n, p)
    results[x + 1] <- results[x + 1] + 1
  }
  results / iter
}


#' Simulate a Hypergeometric Distribution
#'
#' This function simulates sampling from a finite population without replacement
#' and calculates the relative frequencies of the number of successes in each sample.
#'
#' @param iter Number of simulations to run
#' @param n Sample size
#' @param N Population size
#' @param r Number of successes in the population
#'
#' @return A vector of relative frequencies for 0 to n successes
#' @export
#'
#' @examples
#' myhyper(iter = 1000, n = 5, N = 20, r = 12)
myhyper <- function(iter = 100, n, N, r) {
  results <- numeric(n + 1)
  for (i in 1:iter) {
    x <- rhyper(1, m = r, n = N - r, k = n)
    results[x + 1] <- results[x + 1] + 1
  }
  results / iter
}


#' Random Sampling Bar Plots
#'
#' This function repeatedly samples numbers from 1 to 10 (with replacement)
#' and creates a bar plot of the relative frequencies for each sample.
#'
#' @param n Number of draws per sample
#' @param iter Number of iterations (plots) to produce
#' @param time Pause in seconds between plots (default 0.5)
#'
#' @return NULL (produces bar plots only)
#' @export
#'
#' @examples
#' mysample(n = 1000, iter = 1, time = 0)
mysample <- function(n, iter = 10, time = 0.5) {
  for (i in 1:iter) {
    s <- sample(1:10, n, replace = TRUE)
    sf <- factor(s, levels = 1:10)
    barplot(
      table(sf) / n,
      beside = TRUE,
      col = rainbow(10),
      main = paste("Example mysample() iteration", i, "n =", n),
      ylim = c(0, 0.2)
    )
    Sys.sleep(time)
  }
}
