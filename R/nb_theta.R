#' Score and second derivative of the NB log-likelihood in theta
#'
#' Implements the formulas in Section 3.1(b) of the project description,
#' with beta (and therefore mu) held fixed:
#'   U(theta) = sum a_i [ digamma(y_i + theta) - digamma(theta) + log(theta)
#'                        + 1 - log(theta + mu_i)
#'                        - (y_i + theta) / (theta + mu_i) ]
#'   H(theta) = sum a_i [ trigamma(y_i + theta) - trigamma(theta) + 1 / theta
#'                        - 2 / (theta + mu_i)
#'                        + (y_i + theta) / (theta + mu_i)^2 ]
#'
#' @param theta Dispersion parameter (a single number > 0).
#' @param y Vector of counts.
#' @param mu Vector of fitted means.
#' @param weights Prior weights a_i (default: all 1).
#' @return A single number.
#' @noRd
nb_theta_score <- function(theta, y, mu, weights = rep(1, length(y))) {
  sum(weights * (digamma(y + theta) - digamma(theta) + log(theta) + 1 -
                   log(theta + mu) - (y + theta) / (theta + mu)))
}

#' @rdname nb_theta_score
#' @noRd
nb_theta_hess <- function(theta, y, mu, weights = rep(1, length(y))) {
  sum(weights * (trigamma(y + theta) - trigamma(theta) + 1 / theta -
                   2 / (theta + mu) + (y + theta) / (theta + mu)^2))
}

#' Newton-Raphson for theta, with beta held fixed
#'
#' Works on phi = log(theta) so that theta stays positive. By the chain rule
#'   U_phi = theta * U(theta)
#'   H_phi = theta^2 * H(theta) + theta * U(theta)
#' and the update is phi <- phi - U_phi / H_phi.
#'
#' @param y Vector of counts.
#' @param mu Vector of fitted means, held fixed.
#' @param weights Prior weights a_i (default: all 1).
#' @param start Optional starting value for theta. If NULL, a moment
#'   estimate is used.
#' @param tol Convergence tolerance for the change in log(theta).
#' @param maxit Maximum number of Newton-Raphson iterations.
#' @return A list with theta, the number of iterations and a convergence flag.
#' @noRd
nb_theta_newton <- function(y, mu, weights = rep(1, length(y)), start = NULL,
                            tol = 1e-8, maxit = 50) {
  # Moment estimate: Var(y) = mu + mu^2 / theta implies
  # E[(y / mu - 1)^2] is roughly 1 / mu + 1 / theta.
  if (is.null(start)) {
    start <- sum(weights) / sum(weights * (y / mu - 1)^2)
  }
  phi <- log(start)

  converged <- FALSE
  iter <- 0
  while (iter < maxit) {
    iter <- iter + 1
    theta <- exp(phi)

    U <- nb_theta_score(theta, y, mu, weights)
    H <- nb_theta_hess(theta, y, mu, weights)
    U_phi <- theta * U
    H_phi <- theta^2 * H + theta * U

    # Newton-Raphson step. If the curvature has the wrong sign, just move
    # uphill. Never move more than 1 unit on the log scale in one step.
    step <- if (is.finite(H_phi) && H_phi < 0) -U_phi / H_phi else sign(U_phi)
    step <- max(min(step, 1), -1)
    phi <- phi + step

    if (abs(step) <= tol) {
      converged <- TRUE
      break
    }
  }

  list(theta = exp(phi), iter = iter, converged = converged)
}

