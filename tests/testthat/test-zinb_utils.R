test_that("zinb_loglik matches the mixture density written out directly", {
  set.seed(6)
  n <- 80
  mu <- runif(n, 0.5, 6)
  eta_zero <- rnorm(n, -0.5, 1)
  theta <- 1.7
  y <- rnbinom(n, size = theta, mu = mu)
  y[runif(n) < plogis(eta_zero)] <- 0

  p <- plogis(eta_zero)
  dens <- ifelse(y == 0,
                 p + (1 - p) * dnbinom(0, size = theta, mu = mu),
                 (1 - p) * dnbinom(y, size = theta, mu = mu))
  expect_equal(zinb_loglik(y, mu, eta_zero, theta), sum(log(dens)),
               tolerance = 1e-10)
})

test_that("zinb_loglik and zinb_estep agree with pscl::zeroinfl", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  ref <- pscl::zeroinfl(art ~ fem + mar + kid5 + phd + ment,
                        data = bioChemists, dist = "negbin")
  y <- bioChemists$art
  mu <- predict(ref, type = "count")
  p <- predict(ref, type = "zero")

  # Log-likelihood at the zeroinfl estimates.
  expect_equal(zinb_loglik(y, mu, qlogis(p), ref$theta),
               as.numeric(logLik(ref)), tolerance = 1e-8)

  # tau_i = P(structural zero | y_i = 0) = pi_i / P(Y_i = 0).
  tau <- zinb_estep(y, mu, qlogis(p), ref$theta)
  p0 <- predict(ref, type = "prob")[, 1]
  expect_equal(unname(tau[y == 0]), unname(p[y == 0] / p0[y == 0]),
               tolerance = 1e-8)
  expect_true(all(tau[y > 0] == 0))
  expect_true(all(tau >= 0 & tau <= 1))
})

test_that("logit_irls matches glm for a 0/1 response", {
  set.seed(7)
  n <- 400
  Z <- cbind(`(Intercept)` = 1, z1 = rnorm(n), z2 = rbinom(n, 1, 0.4))
  yb <- rbinom(n, 1, plogis(-0.3 + 0.8 * Z[, 2] - 0.5 * Z[, 3]))

  fit <- logit_irls(Z, yb)
  ref <- glm(yb ~ Z - 1, family = binomial,
             control = glm.control(epsilon = 1e-12, maxit = 100))
  expect_true(fit$converged)
  expect_equal(unname(fit$coefficients), unname(coef(ref)), tolerance = 1e-7)
})

test_that("logit_irls matches glm for a fractional response", {
  set.seed(8)
  n <- 400
  Z <- cbind(`(Intercept)` = 1, z1 = rnorm(n))
  tau <- plogis(0.2 + 0.9 * Z[, 2] + rnorm(n, 0, 0.5))

  # quasibinomial has the same estimating equations and accepts
  # non-integer responses without warnings.
  fit <- logit_irls(Z, tau)
  ref <- glm(tau ~ Z - 1, family = quasibinomial,
             control = glm.control(epsilon = 1e-12, maxit = 100))
  expect_equal(unname(fit$coefficients), unname(coef(ref)), tolerance = 1e-7)
})
