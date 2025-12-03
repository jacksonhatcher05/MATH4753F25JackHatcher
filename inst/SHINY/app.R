# Shiny MLE Explorer: Demonstrates MLE for >=5 univariate distributions
# Save this file as app.R and run with: shiny::runApp('.')

library(shiny)
library(ggplot2)
library(gridExtra)
library(stats)

# Negative log-likelihood functions ------------------------------------------------
negloglik_normal <- function(par, x) {
  mu <- par[1]
  sigma <- abs(par[2])
  n <- length(x)
  -sum(dnorm(x, mean = mu, sd = sigma, log = TRUE))
}

negloglik_exponential <- function(par, x) {
  rate <- abs(par[1])
  -sum(dexp(x, rate = rate, log = TRUE))
}

negloglik_poisson <- function(par, x) {
  lambda <- abs(par[1])
  -sum(dpois(x, lambda = lambda, log = TRUE))
}

negloglik_binomial <- function(par, counts, trials) {
  p <- plogis(par[1]) # map to (0,1) using logistic
  -sum(dbinom(counts, size = trials, prob = p, log = TRUE))
}

negloglik_gamma <- function(par, x) {
  shape <- abs(par[1])
  rate <- abs(par[2])
  -sum(dgamma(x, shape = shape, rate = rate, log = TRUE))
}

# Helper for MLE with optim and approximate SE from Hessian -------------------------
mle_optim <- function(negloglik_fn, start_par, lower = NULL, upper = NULL, method = "BFGS", hessian = TRUE, ...) {
  opt <- optim(start_par, negloglik_fn, method = method, hessian = hessian, control = list(maxit = 10000), ...)
  est <- opt$par
  # map transforms for binomial p
  if (!is.null(opt$convergence) && opt$convergence != 0) warning("optim did not converge (code: ", opt$convergence, ")")
  se <- NA
  if (!is.null(opt$hessian) && all(!is.na(opt$hessian))) {
    invh <- try(solve(opt$hessian), silent = TRUE)
    if (!inherits(invh, "try-error")) se <- sqrt(diag(invh))
  }
  list(opt = opt, est = est, se = se)
}

# UI -------------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Shiny MLE Explorer — Univariate Distributions"),
  sidebarLayout(
    sidebarPanel(
      selectInput("dist", "Choose distribution:",
                  choices = c("Normal", "Exponential", "Poisson", "Binomial", "Gamma"), selected = "Normal"),
      uiOutput("param_ui"),
      numericInput("n", "Sample size (or number of observations)", value = 100, min = 1, step = 1),
      numericInput("seed", "Random seed (0 = random)", value = 0, step = 1),
      actionButton("generate", "Generate data & compute MLE"),
      hr(),
      checkboxInput("show_ll", "Show log-likelihood surface", value = TRUE),
      checkboxInput("show_profile", "Show profile likelihood (1D) when available", value = TRUE),
      hr(),
      helpText("This app generates data from a selected univariate distribution, computes MLEs (using analytical formulas where available or numerical optimization), and displays plots and diagnostics.")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Summary & Estimates", verbatimTextOutput("summary")),
        tabPanel("Data & Fits", plotOutput("dataPlot", height = "480px")),
        tabPanel("Log-likelihood", plotOutput("llPlot", height = "480px")),
        tabPanel("Profile", plotOutput("profilePlot", height = "480px")),
        tabPanel("Code & Notes",
                 h4("Notes"),
                 p("This single-file app shows MLE for: Normal (mu,sigma), Exponential (rate), Poisson (lambda), Binomial (p), Gamma (shape,rate).\nMLEs for Poisson, Binomial, Exponential, and Normal have closed forms; gamma requires numerical optimization.")
        )
      )
    )
  )
)

# Server ---------------------------------------------------------------------------
server <- function(input, output, session) {
  # Dynamic UI for parameters depending on distribution
  output$param_ui <- renderUI({
    switch(input$dist,
           "Normal" = tagList(
             numericInput("mu", "True mu:", value = 0),
             numericInput("sigma", "True sigma (>0):", value = 1, min = 1e-6)
           ),
           "Exponential" = numericInput("rate", "True rate (lambda) (>0):", value = 1, min = 1e-6),
           "Poisson" = numericInput("lambda", "True lambda (>0):", value = 3, min = 0),
           "Binomial" = tagList(
             numericInput("p", "True p (0-1):", value = 0.3, min = 0, max = 1, step = 0.01),
             numericInput("trials", "Trials per observation (n for each binomial obs):", value = 10, min = 1, step = 1)
           ),
           "Gamma" = tagList(
             numericInput("shape", "True shape (>0):", value = 2, min = 1e-6),
             numericInput("rate_g", "True rate (>0):", value = 1, min = 1e-6)
           )
    )
  })

  # Reactive: generated data and MLE results
  sim <- eventReactive(input$generate, {
    if (input$seed != 0) set.seed(input$seed)
    n <- input$n
    dist <- input$dist
    res <- list()
    if (dist == "Normal") {
      x <- rnorm(n, mean = input$mu, sd = input$sigma)
      # closed form MLEs
      mu_hat <- mean(x)
      sigma_hat <- sqrt(sum((x - mu_hat)^2)/n)
      # approximate SEs
      se_mu <- sigma_hat / sqrt(n)
      se_sigma <- sigma_hat / sqrt(2*n)
      res <- list(x = x, mle = c(mu = mu_hat, sigma = sigma_hat), se = c(mu = se_mu, sigma = se_sigma), method = "closed form")
    } else if (dist == "Exponential") {
      x <- rexp(n, rate = input$rate)
      rate_hat <- 1/mean(x)
      se_rate <- rate_hat / sqrt(n)
      res <- list(x = x, mle = c(rate = rate_hat), se = c(rate = se_rate), method = "closed form")
    } else if (dist == "Poisson") {
      x <- rpois(n, lambda = input$lambda)
      lambda_hat <- mean(x)
      se_lambda <- sqrt(lambda_hat / n)
      res <- list(x = x, mle = c(lambda = lambda_hat), se = c(lambda = se_lambda), method = "closed form")
    } else if (dist == "Binomial") {
      trials <- input$trials
      # simulate binomial observations (each obs is count of successes in 'trials')
      counts <- rbinom(n, size = trials, prob = input$p)
      p_hat <- mean(counts)/trials
      se_p <- sqrt(p_hat*(1-p_hat)/(n*trials))
      res <- list(counts = counts, trials = trials, mle = c(p = p_hat), se = c(p = se_p), method = "closed form")
    } else if (dist == "Gamma") {
      x <- rgamma(n, shape = input$shape, rate = input$rate_g)
      # initial guesses via method of moments
      m <- mean(x); v <- var(x); shape0 <- m^2 / v; rate0 <- m / v
      neglog <- function(par) negloglik_gamma(par, x)
      opt <- mle_optim(neglog, start_par = c(shape0, rate0), method = "BFGS")
      est <- abs(opt$est)
      se <- opt$se
      res <- list(x = x, mle = c(shape = est[1], rate = est[2]), se = se, method = "optim", opt = opt$opt)
    }
    res
  })

  output$summary <- renderPrint({
    req(sim())
    dist <- input$dist
    cat("Distribution:", dist, "\n")
    if (dist == "Binomial") {
      cat("True p:", input$p, " trials per obs:", input$trials, "\n")
      cat("MLE (p):", round(sim()$mle["p"], 5), "\n")
      cat("Approx SE:", round(sim()$se["p"], 5), "\n")
      cat("Method: closed form (sample mean / trials)\n")
      cat("Sample of counts (first 10):", paste(head(sim()$counts, 10), collapse = ", "))
    } else {
      print(sim())
    }
  })

  output$dataPlot <- renderPlot({
    req(sim())
    dist <- input$dist
    if (dist == "Binomial") {
      df <- data.frame(counts = sim()$counts)
      ggplot(df, aes(x = counts)) +
        geom_bar(aes(y = ..prop..), stat = "count") +
        geom_vline(xintercept = sim()$mle["p"]*sim()$trials, linetype = "dashed") +
        labs(title = "Histogram of Binomial counts", y = "Proportion") +
        theme_minimal()
    } else {
      x <- sim()$x
      df <- data.frame(x = x)
      p <- ggplot(df, aes(x = x)) + geom_histogram(aes(y = ..density..), bins = 30, fill = NA) + theme_minimal()
      # overlay fit/density
      if (dist == "Normal") {
        p <- p + stat_function(fun = dnorm, args = list(mean = sim()$mle["mu"], sd = sim()$mle["sigma"])) + labs(title = "Histogram & Normal fit (MLE)")
      } else if (dist == "Exponential") {
        p <- p + stat_function(fun = dexp, args = list(rate = sim()$mle["rate"])) + labs(title = "Histogram & Exponential fit (MLE)")
      } else if (dist == "Poisson") {
        # show empirical barplot of counts
        df2 <- data.frame(x = x)
        return(ggplot(df2, aes(x = factor(x))) + geom_bar(aes(y = ..prop..), stat = "count") + labs(title = "Poisson counts (proportion)") + theme_minimal())
      } else if (dist == "Gamma") {
        p <- p + stat_function(fun = dgamma, args = list(shape = sim()$mle["shape"], rate = sim()$mle["rate"])) + labs(title = "Histogram & Gamma fit (MLE)")
      }
      p
    }
  })

  output$llPlot <- renderPlot({
    req(sim())
    if (!input$show_ll) return(NULL)
    dist <- input$dist
    n <- input$n

    if (dist == "Normal") {
      x <- sim()$x
      # grid for mu and sigma
      mu_seq <- seq(mean(x) - 2*sd(x), mean(x) + 2*sd(x), length.out = 80)
      sigma_seq <- seq(max(sd(x)/5, 1e-3), sd(x)*2, length.out = 80)
      grid <- expand.grid(mu = mu_seq, sigma = sigma_seq)
      grid$ll <- mapply(function(mu, sigma) sum(dnorm(x, mean = mu, sd = sigma, log = TRUE)), grid$mu, grid$sigma)
      ggplot(grid, aes(x = mu, y = sigma, z = ll)) + geom_contour_filled() + geom_point(aes(x = sim()$mle["mu"], y = sim()$mle["sigma"]), color = "red", size = 3) + labs(title = "Log-likelihood surface (Normal)") + theme_minimal()
    } else if (dist == "Gamma") {
      x <- sim()$x
      shape_seq <- seq(max(0.01, sim()$mle["shape"]*0.3), sim()$mle["shape"]*2.5, length.out = 80)
      rate_seq <- seq(max(0.01, sim()$mle["rate"]*0.3), sim()$mle["rate"]*2.5, length.out = 80)
      grid <- expand.grid(shape = shape_seq, rate = rate_seq)
      grid$ll <- mapply(function(sh, ra) sum(dgamma(x, shape = sh, rate = ra, log = TRUE)), grid$shape, grid$rate)
      ggplot(grid, aes(x = shape, y = rate, z = ll)) + geom_contour_filled() + geom_point(aes(x = sim()$mle["shape"], y = sim()$mle["rate"]), color = "red", size = 3) + labs(title = "Log-likelihood surface (Gamma)") + theme_minimal()
    } else if (dist %in% c("Exponential", "Poisson", "Binomial")) {
      # 1D log-likelihood plots
      if (dist == "Exponential") {
        x <- sim()$x
        rate_seq <- seq(max(1e-3, sim()$mle["rate"]*0.3), sim()$mle["rate"]*2.5, length.out = 200)
        ll <- sapply(rate_seq, function(r) sum(dexp(x, rate = r, log = TRUE)))
        df <- data.frame(rate = rate_seq, ll = ll)
        ggplot(df, aes(x = rate, y = ll)) + geom_line() + geom_vline(xintercept = sim()$mle["rate"], linetype = "dashed") + labs(title = "Log-likelihood (Exponential)") + theme_minimal()
      } else if (dist == "Poisson") {
        x <- sim()$x
        lambda_seq <- seq(max(0, sim()$mle["lambda"]*0.3), sim()$mle["lambda"]*2.5 + 1, length.out = 200)
        ll <- sapply(lambda_seq, function(l) sum(dpois(x, lambda = l, log = TRUE)))
        df <- data.frame(lambda = lambda_seq, ll = ll)
        ggplot(df, aes(x = lambda, y = ll)) + geom_line() + geom_vline(xintercept = sim()$mle["lambda"], linetype = "dashed") + labs(title = "Log-likelihood (Poisson)") + theme_minimal()
      } else if (dist == "Binomial") {
        counts <- sim()$counts
        trials <- sim()$trials
        p_seq <- seq(0.001, 0.999, length.out = 200)
        ll <- sapply(p_seq, function(p) sum(dbinom(counts, size = trials, prob = p, log = TRUE)))
        df <- data.frame(p = p_seq, ll = ll)
        ggplot(df, aes(x = p, y = ll)) + geom_line() + geom_vline(xintercept = sim()$mle["p"], linetype = "dashed") + labs(title = "Log-likelihood (Binomial)") + theme_minimal()
      }
    }
  })

  output$profilePlot <- renderPlot({
    req(sim())
    if (!input$show_profile) return(NULL)
    dist <- input$dist
    if (dist == "Normal") {
      x <- sim()$x
      mu_seq <- seq(mean(x) - 2*sd(x), mean(x) + 2*sd(x), length.out = 200)
      prof <- sapply(mu_seq, function(mu) {
        sigma_hat <- sqrt(sum((x - mu)^2)/length(x))
        sum(dnorm(x, mean = mu, sd = sigma_hat, log = TRUE))
      })
      df <- data.frame(mu = mu_seq, prof = prof)
      ggplot(df, aes(x = mu, y = prof)) + geom_line() + geom_vline(xintercept = sim()$mle["mu"], linetype = "dashed") + labs(title = "Profile log-likelihood for mu (Normal)") + theme_minimal()
    } else if (dist == "Gamma") {
      x <- sim()$x
      shape_seq <- seq(max(0.01, sim()$mle["shape"]*0.3), sim()$mle["shape"]*2.5, length.out = 200)
      prof <- sapply(shape_seq, function(sh) {
        # maximize over rate for each shape
        rates <- optimize(function(r) -sum(dgamma(x, shape = sh, rate = r, log = TRUE)), interval = c(1e-6, sim()$mle["rate"]*5))
        val <- sum(dgamma(x, shape = sh, rate = rates$minimum, log = TRUE))
        val
      })
      df <- data.frame(shape = shape_seq, prof = prof)
      ggplot(df, aes(x = shape, y = prof)) + geom_line() + geom_vline(xintercept = sim()$mle["shape"], linetype = "dashed") + labs(title = "Profile log-likelihood for shape (Gamma)") + theme_minimal()
    } else if (dist == "Exponential") {
      x <- sim()$x
      rate_seq <- seq(max(1e-3, sim()$mle["rate"]*0.3), sim()$mle["rate"]*2.5, length.out = 200)
      prof <- sapply(rate_seq, function(r) sum(dexp(x, rate = r, log = TRUE)))
      df <- data.frame(rate = rate_seq, prof = prof)
      ggplot(df, aes(x = rate, y = prof)) + geom_line() + geom_vline(xintercept = sim()$mle["rate"], linetype = "dashed") + labs(title = "Profile log-likelihood (Exponential)") + theme_minimal()
    } else if (dist == "Poisson") {
      x <- sim()$x
      lambda_seq <- seq(max(0.001, sim()$mle["lambda"]*0.3), sim()$mle["lambda"]*2.5 + 1, length.out = 200)
      prof <- sapply(lambda_seq, function(l) sum(dpois(x, lambda = l, log = TRUE)))
      df <- data.frame(lambda = lambda_seq, prof = prof)
      ggplot(df, aes(x = lambda, y = prof)) + geom_line() + geom_vline(xintercept = sim()$mle["lambda"], linetype = "dashed") + labs(title = "Profile log-likelihood (Poisson)") + theme_minimal()
    } else if (dist == "Binomial") {
      counts <- sim()$counts
      trials <- sim()$trials
      p_seq <- seq(0.001, 0.999, length.out = 200)
      prof <- sapply(p_seq, function(p) sum(dbinom(counts, size = trials, prob = p, log = TRUE)))
      df <- data.frame(p = p_seq, prof = prof)
      ggplot(df, aes(x = p, y = prof)) + geom_line() + geom_vline(xintercept = sim()$mle["p"], linetype = "dashed") + labs(title = "Profile log-likelihood (Binomial)") + theme_minimal()
    }
  })

}

# Run the app
shinyApp(ui, server)
