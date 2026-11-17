test_that("theta score and second derivative match numerical derivatives", {
  set.seed(4)
  mu <- runif(60, 0.5, 6)
  y <- rnbinom(60, size = 1.8, mu = mu)
  a <- runif(60, 0.3, 1)
  theta <- 2.1
  h <- 1e-5

  num_score <- (nb_loglik(y, mu, theta + h, a) -
                  nb_loglik(y, mu, theta - h, a)) / (2 * h)
  expect_equal(nb_theta_score(theta, y, mu, a), num_score, tolerance = 1e-6)

  num_hess <- (nb_theta_score(theta + h, y, mu, a) -
                 nb_theta_score(theta - h, y, mu, a)) / (2 * h)
  expect_equal(nb_theta_hess(theta, y, mu, a), num_hess, tolerance = 1e-6)
})

test_that("nb_theta_newton reproduces theta and its SE from glm.nb", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  ref <- MASS::glm.nb(art ~ fem + mar + kid5 + phd + ment, data = bioChemists,
                      control = glm.control(epsilon = 1e-12, maxit = 100))
  y <- bioChemists$art
  mu <- fitted(ref)

  fit <- nb_theta_newton(y, mu)
  expect_true(fit$converged)
  expect_equal(fit$theta, ref$theta, tolerance = 1e-5)

  # The score is zero at the maximum.
  expect_lt(abs(nb_theta_score(fit$theta, y, mu)), 1e-6)

  # SE(theta) = sqrt(-1 / H(theta)).
  se <- sqrt(-1 / nb_theta_hess(fit$theta, y, mu))
  expect_equal(se, ref$SE.theta, tolerance = 1e-4)
})

