#' ZINB log-likelihood as a function of one parameter vector
#'
#' Wraps zinb_loglik() so that it takes all parameters in a single vector
#' par = (beta, gamma, log(theta)), the form needed to differentiate it.
#'
#' @param par Parameter vector (beta, gamma, log(theta)).
#' @param X Count-model design matrix.
#' @param Z Zero-model design matrix.
#' @param y Vector of counts.
#' @param offset Count-model offset.
#' @return A single number: the observed-data log-likelihood.
#' @noRd
zinb_loglik_par <- function(par, X, Z, y, offset) {
  p <- ncol(X)
  q <- ncol(Z)
  beta <- par[seq_len(p)]
  gamma <- par[p + seq_len(q)]
  theta <- exp(par[p + q + 1])
  mu <- exp(drop(X %*% beta) + offset)
  zinb_loglik(y, mu, drop(Z %*% gamma), theta)
}

#' Variance-covariance matrix of the ZINB estimates
#'
#' EM does not produce standard errors. Following Section 3.2 of the
#' project description, we compute the observed information matrix of the
#' observed-data log-likelihood with respect to (beta, gamma, log(theta))
#' at the MLE and invert it:
#'   Var(estimates) = [ -Hessian of log-likelihood ]^{-1}.
#' The Hessian is either derived analytically (zinb_hessian(), the default)
#' or obtained numerically with numDeriv::hessian(), which differentiates
#' our own log-likelihood (it does not optimize anything).
#'
#' @param beta,gamma,theta Estimates at convergence.
#' @param hessian "analytic" or "numerical".
#' @inheritParams zinb_loglik_par
#' @return The variance-covariance matrix, with rows and columns named
#'   count_*, zero_* and Log(theta), as in pscl::zeroinfl().
#' @noRd
zinb_vcov <- function(beta, gamma, theta, X, Z, y, offset,
                      hessian = c("analytic", "numerical")) {
  hessian <- match.arg(hessian)
  par <- c(beta, gamma, log(theta))
  H <- if (hessian == "analytic") {
    zinb_hessian(par, X, Z, y, offset)
  } else {
    numDeriv::hessian(zinb_loglik_par, par, X = X, Z = Z, y = y,
                      offset = offset)
  }
  V <- solve(-H)
  nms <- c(paste0("count_", colnames(X)), paste0("zero_", colnames(Z)),
           "Log(theta)")
  dimnames(V) <- list(nms, nms)
  V
}
