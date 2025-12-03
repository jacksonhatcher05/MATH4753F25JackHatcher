#' Find number of tickets to sell (discrete and normal approx)
#'
#' \code{ntickets} computes the number of tickets to sell given
#' the number of seats \code{N}, show probability \code{p}, and target
#' overbooking probability \code{gamma}. It returns both a discrete
#' (binomial) solution and a normal-approximation solution, and plots
#' objective functions for both.
#'
#' @param N Integer. Number of seats on the flight.
#' @param gamma Numeric in (0,1). Acceptable probability that the plane
#'   is truly overbooked (more people show than seats).
#' @param p Numeric in (0,1). Probability that an individual ticket-holder shows up.
#' @return A named list with elements:
#'   \item{nd}{integer: number of tickets from exact discrete/binomial calculation}
#'   \item{nc}{integer: number of tickets from normal approximation (continuity corrected)}
#'   \item{N}{the seats input}
#'   \item{p}{show probability input}
#'   \item{gamma}{target overbooking probability input}
#' @examples
#' ntickets(N = 400, gamma = 0.02, p = 0.95)
#' @export
ntickets <- function(N, gamma, p) {
  if (!is.numeric(N) || length(N) != 1 || N <= 0) stop("N must be a positive integer (or numeric).")
  if (!is.numeric(p) || p <= 0 || p >= 1) stop("p must be in (0,1).")
  if (!is.numeric(gamma) || gamma <= 0 || gamma >= 1) stop("gamma must be in (0,1).")

  N <- as.integer(N)

  ## 1. Discrete solution (exact binomial)
  # find smallest integer n >= N such that P(#shows > N) <= gamma
  # i.e., 1 - pbinom(N, size = n, prob = p) <= gamma
  f_disc <- function(n) {
    1 - gamma - pbinom(N, size = n, prob = p)
  }

  # search upward from N until condition met; expand upper bound if needed
  max_increment <- 10000
  upper <- N
  found_nd <- NA_integer_
  while (upper <= N + max_increment) {
    val <- f_disc(upper)
    if (val <= 0) { found_nd <- upper; break }
    upper <- upper + 1L
  }
  if (is.na(found_nd)) {
    stop("Discrete search failed to find nd within reasonable search range. Try increasing search bound.")
  }
  nd <- found_nd

  ## 2. Continuous solution using normal approx (with continuity correction)
  # For Bin(n,p) approx Normal(mean = n*p, sd = sqrt(n*p*(1-p)))
  # We solve 1 - gamma - Phi( (N + 0.5 - n*p) / sqrt(n*p*(1-p)) ) = 0
  f_cont <- function(n) {
    # allow n as real (uniroot operates on reals)
    # but for sqrt argument we must have n*p*(1-p) > 0
    mu <- n * p
    sigma2 <- n * p * (1 - p)
    if (sigma2 <= 0) return(1.0) # large positive so root not found here
    z <- (N + 0.5 - mu) / sqrt(sigma2)
    1 - gamma - pnorm(z)
  }

  # Find bracket for uniroot: function f_cont(n) decreases in n generally; evaluate at N and a high upper
  left <- max(N, 1)
  right <- left + 2000
  # expand right until sign change or too large
  fl <- f_cont(left)
  fr <- f_cont(right)
  attempts <- 0
  while (!(fl * fr <= 0) && attempts < 20) {
    right <- right + 2000
    fr <- f_cont(right)
    attempts <- attempts + 1
  }
  if (fl * fr > 0) {
    # fallback: search by minimizing absolute value over a candidate grid and pick which.min
    cand_ns <- seq(from = N, to = N + 5000, by = 1)
    obj_vals <- vapply(cand_ns, function(n) abs(f_cont(n)), numeric(1))
    n_cont_real <- cand_ns[which.min(obj_vals)]
  } else {
    # uniroot to find real-valued root
    un <- uniroot(f_cont, lower = left, upper = right, tol = .Machine$double.eps^0.5, maxiter = 1000)
    n_cont_real <- un$root
  }

  # round up to next integer ticket count (conservative)
  nc <- ceiling(n_cont_real)

  ## 3. Plot objective functions vs n
  # reasonable plotting window: from N to max(N + 50, nd + 50, nc + 50)
  upper_plot <- max(N + 50, nd + 50, nc + 50)
  ns <- seq(N, upper_plot, by = 1)

  # discrete objective: g_disc(n) = 1 - gamma - pbinom(N, size=n, prob=p)
  g_disc_vals <- vapply(ns, function(n) 1 - gamma - pbinom(N, size = n, prob = p), numeric(1))

  # continuous objective: g_cont(n) = 1 - gamma - pnorm((N + 0.5 - n*p)/sqrt(n*p*(1-p)))
  g_cont_vals <- vapply(ns, function(n) {
    sd2 <- n * p * (1 - p)
    if (sd2 <= 0) return(NA_real_)
    z <- (N + 0.5 - n * p) / sqrt(sd2)
    1 - gamma - pnorm(z)
  }, numeric(1))

  # Plot discrete objective
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  # open a 2-panel plot window
  graphics::par(mfrow = c(2, 1), mar = c(4.2, 4.2, 2, 1))
  plot(ns, g_disc_vals,
       type = "b", pch = 19, cex = 0.7,
       xlab = "n (tickets sold)", ylab = "objective (1 - gamma - P(shows <= N))",
       main = paste0("Discrete objective (Binomial)  — N=", N, ", p=", p, ", gamma=", gamma))
  abline(h = 0, col = "gray40", lty = 2)
  # mark nd
  if (nd %in% ns) points(nd, g_disc_vals[which(ns == nd)], pch = 21, bg = "red", cex = 2)
  legend("topright", legend = paste0("nd = ", nd), bty = "n")

  # Plot continuous objective
  plot(ns, g_cont_vals,
       type = "b", pch = 19, cex = 0.7,
       xlab = "n (tickets sold)", ylab = "objective (1 - gamma - Phi(...))",
       main = "Continuous objective (Normal approx, continuity correction)")
  abline(h = 0, col = "gray40", lty = 2)
  if (nc %in% ns) {
    idx_nc <- which(ns == nc)
    points(nc, g_cont_vals[idx_nc], pch = 21, bg = "blue", cex = 2)
  } else {
    # mark approximate real root if not integer in ns
    points(nc, NA) # nothing but keep safe
  }
  legend("topright", legend = paste0("nc = ", nc), bty = "n")

  # restore par
  par(old_par)

  result <- list(nd = nd, nc = nc, N = N, p = p, gamma = gamma)
  class(result) <- "ntickets_result"

  # print result (so it is visible in knit ted HTML)
  print(result)
  invisible(result)
}

