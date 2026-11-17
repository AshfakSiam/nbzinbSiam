# Simulation study: does zinbreg() recover known parameters?
#
# Following Section 5.1 of the project description, we simulate R = 100
# data sets of size n = 500 from a ZINB model with one covariate in each
# part and known (beta, gamma, theta), fit each with zinbreg(), and report
#   * the bias of each estimate,
#   * the mean of the estimated standard errors vs. the Monte Carlo SD
#     (the standard deviation of the estimates across data sets),
#   * the empirical coverage of the 95% Wald intervals.
#
# Run from the package's top-level folder, after installing the package:
#   source("inst/benchmark/simulation.R")
# The results are saved to inst/benchmark/simulation_results.rds, which the
# report reads.

library(nbzinbSiam)

set.seed(2026)
R <- 100
n <- 500

# True parameters.
beta <- c(`(Intercept)` = 1.0, x = 0.5)    # count model: log mu = 1 + 0.5 x
gamma <- c(`(Intercept)` = -1.0, z = 0.8)  # zero model: logit pi = -1 + 0.8 z
theta <- 2
truth <- c(beta, gamma, log(theta))
names(truth) <- c("count_(Intercept)", "count_x",
                  "zero_(Intercept)", "zero_z", "Log(theta)")

simulate_zinb <- function(n) {
  x <- rnorm(n)
  z <- rnorm(n)
  mu <- exp(beta[1] + beta[2] * x)
  p <- plogis(gamma[1] + gamma[2] * z)
  structural_zero <- runif(n) < p
  y <- ifelse(structural_zero, 0, rnbinom(n, size = theta, mu = mu))
  data.frame(y = y, x = x, z = z)
}

est <- se <- matrix(NA_real_, R, length(truth),
                    dimnames = list(NULL, names(truth)))
converged <- logical(R)
iterations <- integer(R)
zero_share <- numeric(R)

time_total <- system.time({
  for (r in seq_len(R)) {
    d <- simulate_zinb(n)
    zero_share[r] <- mean(d$y == 0)
    fit <- tryCatch(zinbreg(y ~ x | z, data = d), error = function(e) NULL)
    if (is.null(fit)) next
    est[r, ] <- c(coef(fit), log(fit$theta))
    se[r, ] <- c(sqrt(diag(vcov(fit))), fit$SE.logtheta)
    converged[r] <- fit$converged
    iterations[r] <- fit$iterations
  }
})[["elapsed"]]

ok <- converged & stats::complete.cases(est)
lower <- est - qnorm(0.975) * se
upper <- est + qnorm(0.975) * se
covered <- sweep(lower, 2, truth, `<=`) & sweep(upper, 2, truth, `>=`)

summary_table <- data.frame(
  Parameter = names(truth),
  True = truth,
  Mean = colMeans(est[ok, , drop = FALSE]),
  Bias = colMeans(est[ok, , drop = FALSE]) - truth,
  MC_SD = apply(est[ok, , drop = FALSE], 2, sd),
  Mean_SE = colMeans(se[ok, , drop = FALSE]),
  Coverage = colMeans(covered[ok, , drop = FALSE]),
  row.names = NULL, check.names = FALSE
)

simulation <- list(
  summary = summary_table,
  settings = list(R = R, n = n, beta = beta, gamma = gamma, theta = theta,
                  seed = 2026),
  estimates = est, std_errors = se,
  converged = converged, iterations = iterations, zero_share = zero_share,
  time_total = time_total,
  session = sessionInfo()
)

out_file <- file.path("inst", "benchmark", "simulation_results.rds")
saveRDS(simulation, out_file)

print(summary_table, digits = 3)
cat("\nConverged:", sum(converged), "of", R, "data sets",
    "\nMean EM iterations:", round(mean(iterations[ok]), 1),
    "\nMean share of zeros:", round(mean(zero_share), 3),
    "\nTotal time:", round(time_total, 1), "seconds",
    "\nSaved to", out_file, "\n")
