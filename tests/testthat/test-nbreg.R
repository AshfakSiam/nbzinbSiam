test_that("nbreg reproduces MASS::glm.nb on bioChemists", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  fm <- art ~ fem + mar + kid5 + phd + ment

  fit <- nbreg(fm, data = bioChemists)
  ref <- MASS::glm.nb(fm, data = bioChemists,
                      control = glm.control(epsilon = 1e-10, maxit = 100))

  expect_s3_class(fit, "nbreg")
  expect_true(fit$converged)

  # Tolerances from Table 4 of the project description.
  expect_lt(max(abs(fit$coefficients - coef(ref))), 1e-5)
  expect_lt(max(abs(sqrt(diag(fit$vcov)) - sqrt(diag(vcov(ref))))), 1e-5)
  expect_lt(abs(fit$theta - ref$theta), 1e-4)
  expect_lt(abs(fit$loglik - as.numeric(logLik(ref))), 1e-4)
  expect_lt(abs(fit$SE.theta - ref$SE.theta), 1e-4)
})

test_that("nbreg handles an offset given either way", {
  skip_if_not_installed("MASS")
  set.seed(3)
  n <- 300
  d <- data.frame(x = rnorm(n), expo = runif(n, 0.5, 2))
  d$y <- rnbinom(n, size = 1.5, mu = d$expo * exp(0.3 + 0.6 * d$x))

  fit1 <- nbreg(y ~ x + offset(log(expo)), data = d)
  fit2 <- nbreg(y ~ x, data = d, offset = log(expo))
  ref <- MASS::glm.nb(y ~ x + offset(log(expo)), data = d,
                      control = glm.control(epsilon = 1e-10, maxit = 100))

  expect_lt(max(abs(fit1$coefficients - coef(ref))), 1e-5)
  expect_lt(max(abs(fit2$coefficients - coef(ref))), 1e-5)
  expect_lt(abs(fit1$theta - ref$theta), 1e-4)
})

test_that("nbreg works with a single covariate and with intercept only", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")

  fit1 <- nbreg(art ~ ment, data = bioChemists)
  ref1 <- MASS::glm.nb(art ~ ment, data = bioChemists,
                       control = glm.control(epsilon = 1e-10, maxit = 100))
  expect_lt(max(abs(fit1$coefficients - coef(ref1))), 1e-5)
  expect_lt(abs(fit1$theta - ref1$theta), 1e-4)

  fit0 <- nbreg(art ~ 1, data = bioChemists)
  ref0 <- MASS::glm.nb(art ~ 1, data = bioChemists)
  expect_lt(abs(fit0$coefficients - coef(ref0)), 1e-5)
})

test_that("nbreg gives informative errors for bad input", {
  d <- data.frame(x = 1:6)
  expect_error(nbreg(y ~ x, data = transform(d, y = c(1.5, 2, 3, 4, 5, 6))),
               "non-negative integers")
  expect_error(nbreg(y ~ x, data = transform(d, y = c(-1, 2, 3, 4, 5, 6))),
               "non-negative integers")
  expect_error(nbreg(y ~ x, data = transform(d, y = letters[1:6])),
               "numeric vector of counts")
  expect_error(nbreg(y ~ x, data = transform(d, y = 1:6), control = list()),
               "nbreg_control")
  expect_error(nbreg_control(tol = -1), "positive")
})

test_that("nbreg warns at the Poisson limit and when maxit is reached", {
  set.seed(5)
  d <- data.frame(x = rnorm(400))
  d$y <- rpois(400, exp(1 + 0.3 * d$x))
  expect_warning(fit <- nbreg(y ~ x, data = d), "theta is very large")

  # At the Poisson limit the coefficients equal Poisson regression.
  ref <- glm(y ~ x, data = d, family = poisson)
  expect_lt(max(abs(fit$coefficients - coef(ref))), 1e-5)

  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  expect_warning(nbreg(art ~ fem + ment, data = bioChemists,
                       control = nbreg_control(maxit = 1)),
                 "did not converge")
})

