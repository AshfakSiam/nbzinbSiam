test_that("nb_irls with fixed theta matches glm with negative.binomial family", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  X <- model.matrix(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
  y <- bioChemists$art
  theta <- 2.5

  fit <- nb_irls(X, y, theta)
  ref <- glm(art ~ fem + mar + kid5 + phd + ment, data = bioChemists,
             family = MASS::negative.binomial(theta),
             control = glm.control(epsilon = 1e-14, maxit = 100))

  expect_true(fit$converged)
  expect_equal(fit$coefficients, coef(ref), tolerance = 1e-7)
  expect_equal(fit$loglik,
               sum(dnbinom(y, size = theta, mu = fitted(ref), log = TRUE)),
               tolerance = 1e-9)
})

test_that("nb_irls approaches Poisson regression as theta grows", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  X <- model.matrix(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
  y <- bioChemists$art

  fit <- nb_irls(X, y, theta = 1e7)
  ref <- glm(art ~ fem + mar + kid5 + phd + ment, data = bioChemists,
             family = poisson)

  expect_equal(fit$coefficients, coef(ref), tolerance = 1e-5)
})

test_that("nb_irls handles an offset and prior weights", {
  skip_if_not_installed("MASS")
  set.seed(3)
  n <- 300
  x <- rnorm(n)
  expo <- runif(n, 0.5, 2)
  a <- runif(n, 0.2, 1)
  y <- rnbinom(n, size = 1.5, mu = expo * exp(0.3 + 0.6 * x))
  X <- cbind(`(Intercept)` = 1, x = x)

  fit <- nb_irls(X, y, theta = 1.5, offset = log(expo), weights = a)
  ref <- glm(y ~ x + offset(log(expo)), weights = a,
             family = MASS::negative.binomial(1.5),
             control = glm.control(epsilon = 1e-14, maxit = 100))

  expect_equal(fit$coefficients, coef(ref), tolerance = 1e-7)
})

