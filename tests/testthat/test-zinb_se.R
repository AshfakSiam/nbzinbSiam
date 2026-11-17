test_that("zinbreg standard errors match pscl::zeroinfl", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  fm <- art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment

  fit <- zinbreg(fm, data = bioChemists)
  ref <- pscl::zeroinfl(fm, data = bioChemists, dist = "negbin",
                        control = pscl::zeroinfl.control(reltol = 1e-12))

  # zeroinfl's vcov() covers the coefficients; SE.logtheta is separate.
  se_fit <- sqrt(diag(fit$vcov))
  se_ref <- sqrt(diag(vcov(ref)))
  common <- names(se_ref)
  expect_lt(max(abs(se_fit[common] - se_ref)), 1e-3)
  expect_lt(abs(fit$SE.logtheta - ref$SE.logtheta), 1e-3)
})

test_that("zinbreg vcov is symmetric, positive definite and well named", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  fit <- zinbreg(art ~ fem + ment | ment, data = bioChemists)
  V <- fit$vcov

  expect_equal(V, t(V), tolerance = 1e-8)
  expect_true(all(eigen(V, symmetric = TRUE, only.values = TRUE)$values > 0))
  expect_equal(rownames(V),
               c("count_(Intercept)", "count_femWomen", "count_ment",
                 "zero_(Intercept)", "zero_ment", "Log(theta)"))
  expect_equal(fit$SE.theta, fit$theta * fit$SE.logtheta)
})
