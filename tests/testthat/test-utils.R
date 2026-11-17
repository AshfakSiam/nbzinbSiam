test_that("nb_loglik matches dnbinom", {
  set.seed(1)
  mu <- runif(50, 0.5, 8)
  theta <- 2.3
  y <- rnbinom(50, size = theta, mu = mu)

  expect_equal(nb_loglik(y, mu, theta),
               sum(dnbinom(y, size = theta, mu = mu, log = TRUE)),
               tolerance = 1e-10)

  # with prior weights
  a <- runif(50)
  expect_equal(nb_loglik(y, mu, theta, weights = a),
               sum(a * dnbinom(y, size = theta, mu = mu, log = TRUE)),
               tolerance = 1e-10)
})

test_that("wls_solve matches weighted least squares from lm.wfit", {
  set.seed(2)
  X <- cbind(1, rnorm(40), rnorm(40))
  z <- rnorm(40)
  w <- runif(40, 0.2, 3)

  ref <- unname(stats::lm.wfit(X, z, w)$coefficients)
  expect_equal(unname(wls_solve(X, z, w)), ref, tolerance = 1e-10)
})

