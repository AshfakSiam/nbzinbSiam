# Set up a ZINB problem shared by the tests below.
zinb_test_problem <- function() {
  set.seed(10)
  n <- 300
  d <- data.frame(x = rnorm(n), z = rnorm(n), expo = runif(n, 0.5, 2))
  mu <- d$expo * exp(0.5 + 0.5 * d$x)
  p <- plogis(-0.8 + 0.9 * d$z)
  d$y <- ifelse(runif(n) < p, 0, rnbinom(n, size = 1.5, mu = mu))
  X <- cbind(`(Intercept)` = 1, x = d$x)
  Z <- cbind(`(Intercept)` = 1, z = d$z)
  list(X = X, Z = Z, y = d$y, offset = log(d$expo))
}

test_that("analytic ZINB gradient matches numerical derivatives", {
  pr <- zinb_test_problem()
  # An arbitrary point, not the maximum, so every term is exercised.
  par <- c(0.3, 0.6, -0.5, 0.7, log(1.2))

  num <- numDeriv::grad(zinb_loglik_par, par, X = pr$X, Z = pr$Z,
                        y = pr$y, offset = pr$offset)
  ana <- zinb_score(par, pr$X, pr$Z, pr$y, pr$offset)
  expect_equal(ana, num, tolerance = 1e-6)
})

test_that("analytic ZINB Hessian matches numerical derivatives", {
  pr <- zinb_test_problem()
  par <- c(0.3, 0.6, -0.5, 0.7, log(1.2))

  num <- numDeriv::hessian(zinb_loglik_par, par, X = pr$X, Z = pr$Z,
                           y = pr$y, offset = pr$offset)
  ana <- zinb_hessian(par, pr$X, pr$Z, pr$y, pr$offset)
  expect_equal(unname(ana), num, tolerance = 1e-6)
  expect_equal(ana, t(ana))
})

test_that("the gradient is (nearly) zero at the zinbreg estimates", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  # EM converges slowly, so with the default tolerance the gradient is only
  # about 1e-3. A tighter tolerance brings it close to zero.
  fit <- zinbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists,
                 control = zinbreg_control(tol = 1e-12, maxit = 5000))
  X <- model.matrix(~ fem + mar + kid5 + phd + ment, data = bioChemists)
  par <- c(fit$coefficients$count, fit$coefficients$zero, log(fit$theta))

  g <- zinb_score(par, X, X, bioChemists$art, rep(0, nrow(X)))
  expect_lt(max(abs(g)), 1e-4)
})

test_that("analytic and numerical standard errors agree", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  fm <- art ~ fem + mar + kid5 + phd + ment
  fit_a <- zinbreg(fm, data = bioChemists)
  fit_n <- zinbreg(fm, data = bioChemists,
                   control = zinbreg_control(hessian = "numerical"))

  expect_equal(fit_a$vcov, fit_n$vcov, tolerance = 1e-6)
  expect_equal(fit_a$SE.logtheta, fit_n$SE.logtheta, tolerance = 1e-6)
})
