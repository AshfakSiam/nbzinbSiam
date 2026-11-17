test_that("zinbreg reproduces pscl::zeroinfl on bioChemists", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  fm <- art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment

  fit <- zinbreg(fm, data = bioChemists)
  ref <- pscl::zeroinfl(fm, data = bioChemists, dist = "negbin",
                        control = pscl::zeroinfl.control(reltol = 1e-12))

  expect_s3_class(fit, "zinbreg")
  expect_true(fit$converged)

  # Tolerances from Table 4 of the project description.
  expect_lt(max(abs(fit$coefficients$count - ref$coefficients$count)), 1e-3)
  expect_lt(max(abs(fit$coefficients$zero - ref$coefficients$zero)), 1e-3)
  expect_lt(abs(log(fit$theta) - log(ref$theta)), 1e-3)
  expect_lt(abs(fit$loglik - as.numeric(logLik(ref))), 1e-3)
})

test_that("the EM log-likelihood never decreases", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  fit <- zinbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)

  # Allow only for rounding error in the last digits.
  expect_true(all(diff(fit$loglik_path) > -1e-8))
  expect_equal(length(fit$loglik_path), fit$iterations + 1)
})

test_that("a formula without | uses the same covariates in both parts", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  fit1 <- zinbreg(art ~ fem + ment, data = bioChemists)
  fit2 <- zinbreg(art ~ fem + ment | fem + ment, data = bioChemists)

  expect_equal(fit1$coefficients, fit2$coefficients)
  expect_equal(names(fit1$coefficients$zero),
               c("(Intercept)", "femWomen", "ment"))
})

test_that("zinbreg works with different covariates in the two parts", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  fm <- art ~ fem + kid5 + ment | ment
  fit <- zinbreg(fm, data = bioChemists)
  ref <- pscl::zeroinfl(fm, data = bioChemists, dist = "negbin",
                        control = pscl::zeroinfl.control(reltol = 1e-12))

  expect_lt(max(abs(unlist(fit$coefficients) - unlist(ref$coefficients))),
            1e-3)
})

test_that("zinbreg handles a count-model offset given either way", {
  skip_if_not_installed("pscl")
  set.seed(9)
  n <- 600
  d <- data.frame(x = rnorm(n), z = rnorm(n), expo = runif(n, 0.5, 2))
  mu <- d$expo * exp(0.5 + 0.5 * d$x)
  p <- plogis(-1 + 0.8 * d$z)
  d$y <- ifelse(runif(n) < p, 0, rnbinom(n, size = 2, mu = mu))

  fit1 <- zinbreg(y ~ x + offset(log(expo)) | z, data = d)
  fit2 <- zinbreg(y ~ x | z, data = d, offset = log(expo))
  ref <- pscl::zeroinfl(y ~ x + offset(log(expo)) | z, data = d,
                        dist = "negbin",
                        control = pscl::zeroinfl.control(reltol = 1e-12))

  expect_lt(max(abs(unlist(fit1$coefficients) - unlist(ref$coefficients))),
            1e-3)
  expect_lt(max(abs(unlist(fit2$coefficients) - unlist(ref$coefficients))),
            1e-3)
})

test_that("zinbreg gives informative errors and warnings", {
  d <- data.frame(x = rnorm(10), y = 1:10)
  expect_error(zinbreg(y ~ x, data = d), "no zeros")
  expect_error(zinbreg(y ~ x, data = transform(d, y = c(0, 1.5, 2:9))),
               "non-negative integers")
  expect_error(zinbreg(y ~ x, data = transform(d, y = 0:9), control = list()),
               "zinbreg_control")
  expect_error(zinbreg_control(maxit = 0), "positive")

  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  expect_warning(zinbreg(art ~ ment, data = bioChemists,
                         control = zinbreg_control(maxit = 3)),
                 "did not converge")
})
