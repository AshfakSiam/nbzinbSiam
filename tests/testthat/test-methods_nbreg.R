test_that("nbreg S3 methods agree with glm.nb", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  fm <- art ~ fem + mar + kid5 + phd + ment
  fit <- nbreg(fm, data = bioChemists)
  ref <- MASS::glm.nb(fm, data = bioChemists,
                      control = glm.control(epsilon = 1e-10, maxit = 100))

  expect_equal(coef(fit), coef(ref), tolerance = 1e-6)
  expect_equal(vcov(fit), vcov(ref), tolerance = 1e-6)
  expect_equal(nobs(fit), nobs(ref))

  # logLik carries df, so AIC() and BIC() work.
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_equal(attr(ll, "df"), attr(logLik(ref), "df"))
  expect_equal(AIC(fit), AIC(ref), tolerance = 1e-7)
  expect_equal(BIC(fit), BIC(ref), tolerance = 1e-7)

  # Fitted values and residuals.
  expect_equal(fitted(fit), fitted(ref), tolerance = 1e-6)
  expect_equal(residuals(fit), residuals(ref, type = "response"),
               tolerance = 1e-6)
  expect_equal(residuals(fit, type = "pearson"),
               residuals(ref, type = "pearson"), tolerance = 1e-6)

  # Summary table: estimates, SEs, z values, p-values.
  s <- summary(fit)
  expect_s3_class(s, "summary.nbreg")
  expect_equal(colnames(s$coefficients),
               c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  expect_equal(s$coefficients, summary(ref)$coefficients, tolerance = 1e-5)
})

test_that("predict.nbreg works with and without newdata", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  fm <- art ~ fem + mar + kid5 + phd + ment
  fit <- nbreg(fm, data = bioChemists)
  ref <- MASS::glm.nb(fm, data = bioChemists,
                      control = glm.control(epsilon = 1e-10, maxit = 100))

  expect_equal(predict(fit), fitted(fit))
  expect_equal(predict(fit, type = "link"), log(fitted(fit)))

  # New data, including a factor given with only one of its levels.
  nd <- data.frame(fem = factor("Women", levels = c("Men", "Women")),
                   mar = "Married", kid5 = c(0, 1, 2), phd = 3,
                   ment = c(5, 10, 20))
  expect_equal(unname(predict(fit, newdata = nd)),
               unname(predict(ref, newdata = nd, type = "response")),
               tolerance = 1e-6)
  expect_equal(unname(predict(fit, newdata = nd, type = "link")),
               unname(predict(ref, newdata = nd, type = "link")),
               tolerance = 1e-6)
})

test_that("predict.nbreg uses the offset in newdata", {
  skip_if_not_installed("MASS")
  set.seed(3)
  n <- 300
  d <- data.frame(x = rnorm(n), expo = runif(n, 0.5, 2))
  d$y <- rnbinom(n, size = 1.5, mu = d$expo * exp(0.3 + 0.6 * d$x))
  nd <- data.frame(x = c(-1, 0, 1), expo = c(1, 2, 3))

  fit1 <- nbreg(y ~ x + offset(log(expo)), data = d)
  fit2 <- nbreg(y ~ x, data = d, offset = log(expo))
  expected <- nd$expo * exp(coef(fit1)[1] + coef(fit1)[2] * nd$x)

  expect_equal(unname(predict(fit1, newdata = nd)), unname(expected))
  expect_equal(unname(predict(fit2, newdata = nd)), unname(expected),
               tolerance = 1e-8)
})

test_that("print methods run and return their argument invisibly", {
  set.seed(1)
  d <- data.frame(x = rnorm(100))
  d$y <- rnbinom(100, size = 2, mu = exp(0.5 + 0.8 * d$x))
  fit <- nbreg(y ~ x, data = d)

  expect_output(print(fit), "Theta")
  expect_output(print(summary(fit)), "Pr\\(>\\|z\\|\\)")
  expect_invisible(print(fit))
})

