#' IRLS for the NB regression coefficients, with theta held fixed
#'
#' Implements Section 3.1(a) of the project description. With log link and
#' variance V(mu) = mu + mu^2 / theta, each Fisher-scoring step is a weighted
#' least-squares fit with
#'   working weights   w_i = a_i * theta * mu_i / (theta + mu_i)
#'   working response  z_i = (eta_i - o_i) + (y_i - mu_i) / mu_i
#' where eta_i = x_i' beta + o_i and mu_i = exp(eta_i).
#'
#' @param X Design matrix (n rows, p columns).
#' @param y Vector of counts (length n).
#' @param theta Dispersion parameter, held fixed (a single number > 0).
#' @param offset Offset o_i (default: all 0).
#' @param weights Prior weights a_i (default: all 1).
#' @param start Optional starting value for beta.
#' @param tol Convergence tolerance for the change in beta.
#' @param maxit Maximum number of IRLS iterations.
#' @return A list with coefficients, eta, mu, the final working weights w,
#'   the log-likelihood, the number of iterations and a convergence flag.
#' @noRd
nb_irls <- function(X, y, theta, offset = rep(0, length(y)),
                    weights = rep(1, length(y)), start = NULL,
                    tol = 1e-8, maxit = 50) {
  # Starting value: least squares of log(y + 0.1) on X.
  if (is.null(start)) {
    start <- wls_solve(X, log(y + 0.1) - offset, weights)
  }
  beta <- start
  eta <- drop(X %*% beta) + offset
  mu <- exp(eta)
  ll <- nb_loglik(y, mu, theta, weights)

  converged <- FALSE
  iter <- 0
  while (iter < maxit) {
    iter <- iter + 1

    # Working weights and working response (Section 3.1(a)).
    w <- weights * theta * mu / (theta + mu)
    z <- (eta - offset) + (y - mu) / mu

    # One weighted least-squares step.
    beta_new <- wls_solve(X, z, w)
    eta_new <- drop(X %*% beta_new) + offset
    mu_new <- exp(eta_new)
    ll_new <- nb_loglik(y, mu_new, theta, weights)

    # Step-halving: if the step made the log-likelihood worse (or not a
    # number), move only half as far, up to 10 times.
    half <- 0
    while ((!is.finite(ll_new) || ll_new < ll - 1e-8 * abs(ll)) && half < 10) {
      half <- half + 1
      beta_new <- (beta_new + beta) / 2
      eta_new <- drop(X %*% beta_new) + offset
      mu_new <- exp(eta_new)
      ll_new <- nb_loglik(y, mu_new, theta, weights)
    }

    change <- max(abs(beta_new - beta))
    beta <- beta_new
    eta <- eta_new
    mu <- mu_new
    ll <- ll_new

    if (change <= tol * (1 + max(abs(beta)))) {
      converged <- TRUE
      break
    }
  }

  names(beta) <- colnames(X)
  list(coefficients = beta, eta = eta, mu = mu,
       w = weights * theta * mu / (theta + mu),
       loglik = ll, iter = iter, converged = converged)
}

