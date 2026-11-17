test_that("zinbreg S3 methods agree with zeroinfl", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  fm <- art ~ fem + mar + kid5 + phd + ment
  fit <- zinbreg(fm, data = bioChemists)
  ref <- pscl::zeroinfl(fm, data = bioChemists, dist = "negbin",
                        control = pscl::zeroinfl.control(reltol = 1e-12))

  # Coefficients, named as in zeroinfl.
  expect_equal(names(coef(fit)), names(coef(ref)))
  expect_lt(max(abs(coef(fit) - coef(ref))), 1e-3)
  expect_lt(max(abs(coef(fit, model = "count") -
                      coef(ref, model = "count"))), 1e-3)
  expect_lt(max(abs(coef(fit, model = "zero") -
                      coef(ref, model = "zero"))), 1e-3)

  # vcov covers the coefficients (not Log(theta)), as in zeroinfl.
  expect_equal(dimnames(vcov(fit)), dimnames(vcov(ref)))
  expect_lt(max(abs(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(ref))))), 1e-3)
  expect_equal(dim(vcov(fit, model = "zero")), c(6L, 6L))

  # logLik, AIC, BIC, nobs.
  expect_equal(attr(logLik(fit), "df"), attr(logLik(ref), "df"))
  expect_lt(abs(AIC(fit) - AIC(ref)), 1e-3)
  expect_lt(abs(BIC(fit) - BIC(ref)), 1e-3)
  expect_equal(nobs(fit), ref$n)   # zeroinfl has no nobs() method

  # Fitted values, residuals, predictions.
  expect_lt(max(abs(fitted(fit) - fitted(ref))), 1e-3)
  # Note: residuals() for zeroinfl defaults to Pearson, so ask for "response".
  expect_lt(max(abs(residuals(fit) - residuals(ref, type = "response"))),
            1e-3)
  expect_lt(max(abs(residuals(fit, type = "pearson") -
                      residuals(ref, type = "pearson"))), 1e-3)
  for (tp in c("response", "count", "zero")) {
    expect_lt(max(abs(predict(fit, type = tp) - predict(ref, type = tp))),
              1e-3)
  }
})

test_that("summary.zinbreg has count and zero tables like zeroinfl", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  fm <- art ~ fem + mar + kid5 + phd + ment
  s <- summary(zinbreg(fm, data = bioChemists))
  r <- summary(pscl::zeroinfl(fm, data = bioChemists, dist = "negbin",
                              control = pscl::zeroinfl.control(reltol = 1e-12)))

  expect_s3_class(s, "summary.zinbreg")
  expect_equal(dimnames(s$coefficients$count),
               dimnames(r$coefficients$count))
  expect_equal(dimnames(s$coefficients$zero), dimnames(r$coefficients$zero))
  expect_lt(max(abs(s$coefficients$count[, 1:2] -
                      r$coefficients$count[, 1:2])), 1e-3)
  expect_lt(max(abs(s$coefficients$zero[, 1:2] -
                      r$coefficients$zero[, 1:2])), 1e-3)
})

test_that("predict.zinbreg works with newdata, factors and offsets", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  fm <- art ~ fem + mar + kid5 + phd + ment | fem + ment
  fit <- zinbreg(fm, data = bioChemists)
  ref <- pscl::zeroinfl(fm, data = bioChemists, dist = "negbin",
                        control = pscl::zeroinfl.control(reltol = 1e-12))
  nd <- data.frame(fem = factor("Women", levels = c("Men", "Women")),
                   mar = "Married", kid5 = c(0, 1, 2), phd = 3,
                   ment = c(5, 10, 20))
  for (tp in c("response", "count", "zero")) {
    expect_lt(max(abs(predict(fit, newdata = nd, type = tp) -
                        predict(ref, newdata = nd, type = tp))), 1e-3)
  }
  # Without newdata, predictions are the fitted values.
  expect_equal(predict(fit), fitted(fit))

  # Offset, in the formula or as an argument.
  set.seed(9)
  n <- 600
  d <- data.frame(x = rnorm(n), z = rnorm(n), expo = runif(n, 0.5, 2))
  mu <- d$expo * exp(0.5 + 0.5 * d$x)
  d$y <- ifelse(runif(n) < plogis(-1 + 0.8 * d$z), 0,
                rnbinom(n, size = 2, mu = mu))
  nd2 <- data.frame(x = c(-1, 0, 1), z = 0, expo = c(1, 2, 3))
  f1 <- zinbreg(y ~ x + offset(log(expo)) | z, data = d)
  f2 <- zinbreg(y ~ x | z, data = d, offset = log(expo))
  b <- coef(f1, model = "count")
  expect_equal(unname(predict(f1, newdata = nd2, type = "count")),
               nd2$expo * exp(b[1] + b[2] * nd2$x))
  expect_equal(predict(f1, newdata = nd2), predict(f2, newdata = nd2),
               tolerance = 1e-8)
})

test_that("zinbreg print methods run", {
  set.seed(1)
  n <- 300
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  d$y <- ifelse(runif(n) < plogis(-1 + 0.8 * d$z), 0,
                rnbinom(n, size = 2, mu = exp(1 + 0.5 * d$x)))
  fit <- zinbreg(y ~ x | z, data = d)

  expect_output(print(fit), "Zero-inflation model coefficients")
  expect_output(print(summary(fit)), "Log\\(theta\\)")
  expect_invisible(print(fit))
})
