# Analytic gradient and Hessian of the ZINB observed-data log-likelihood
# with respect to (beta, gamma, phi), where phi = log(theta).
#
# Notation, for observation i:
#   eta_i  = x_i' beta + o_i,  mu_i = exp(eta_i)          (count model)
#   zeta_i = z_i' gamma,       pi_i = plogis(zeta_i)       (zero model)
#   theta  = exp(phi)
#   l_i    = log f_NB(y_i; mu_i, theta)                    (NB log-density)
#
# Derivatives of the NB log-density l_i (any y_i, including 0):
#   d l / d eta        =  theta (y - mu) / (theta + mu)
#   d2 l / d eta2      = -theta mu (theta + y) / (theta + mu)^2
#   d l / d theta      =  U_i: digamma(y + theta) - digamma(theta) + log(theta)
#                         + 1 - log(theta + mu) - (y + theta) / (theta + mu)
#   d2 l / d theta2    =  H_i: trigamma(y + theta) - trigamma(theta)
#                         + 1/theta - 2/(theta + mu) + (y + theta)/(theta + mu)^2
#   d2 l / d eta dtheta =  mu (y - mu) / (theta + mu)^2
# and on the phi = log(theta) scale, by the chain rule:
#   d l / d phi        = theta U_i
#   d2 l / d phi2      = theta^2 H_i + theta U_i
#   d2 l / d eta dphi  = theta mu (y - mu) / (theta + mu)^2
#
# The ZINB contribution is the log of a two-component mixture,
#   ell_i = log( exp(A_i) + exp(B_i) ),
#   A_i = log(pi_i),  B_i = log(1 - pi_i) + l_i,
# with posterior weight tau_i = exp(A_i) / (exp(A_i) + exp(B_i)). For y_i > 0
# the first component is impossible, which is the same as tau_i = 0.
# For such a log-sum the derivatives are
#   grad ell_i = tau grad A + (1 - tau) grad B
#   Hess ell_i = tau Hess A + (1 - tau) Hess B
#                + tau (1 - tau) (grad A - grad B)(grad A - grad B)'.
# With grad A - grad B = (dzeta: 1, deta: -l_eta, dphi: -l_phi), this gives
#   gradient:  beta: (1 - tau) l_eta x,  gamma: (tau - pi) z,
#              phi:  (1 - tau) l_phi
#   Hessian blocks:
#     beta-beta:   [(1 - tau) l_eta,eta + tau (1 - tau) l_eta^2] x x'
#     gamma-gamma: [tau (1 - tau) - pi (1 - pi)] z z'
#     phi-phi:     (1 - tau) l_phi,phi + tau (1 - tau) l_phi^2
#     beta-gamma:  -tau (1 - tau) l_eta  x z'
#     beta-phi:    [(1 - tau) l_eta,phi + tau (1 - tau) l_eta l_phi] x
#     gamma-phi:   -tau (1 - tau) l_phi  z
# Each block is a sum over observations.

#' Per-observation pieces shared by the ZINB gradient and Hessian
#'
#' @inheritParams zinb_loglik_par
#' @return A list with tau, pi and the NB log-density derivatives.
#' @noRd
zinb_deriv_pieces <- function(par, X, Z, y, offset) {
  p <- ncol(X)
  q <- ncol(Z)
  beta <- par[seq_len(p)]
  gamma <- par[p + seq_len(q)]
  theta <- exp(par[p + q + 1])

  eta <- drop(X %*% beta) + offset
  mu <- exp(eta)
  zeta <- drop(Z %*% gamma)
  pi0 <- stats::plogis(zeta)
  tau <- zinb_estep(y, mu, zeta, theta)

  tm <- theta + mu
  l_eta <- theta * (y - mu) / tm
  l_etaeta <- -theta * mu * (theta + y) / tm^2
  U <- digamma(y + theta) - digamma(theta) + log(theta) + 1 - log(tm) -
    (y + theta) / tm
  H <- trigamma(y + theta) - trigamma(theta) + 1 / theta - 2 / tm +
    (y + theta) / tm^2
  l_phi <- theta * U
  l_phiphi <- theta^2 * H + theta * U
  l_etaphi <- theta * mu * (y - mu) / tm^2

  list(tau = tau, pi = pi0, l_eta = l_eta, l_etaeta = l_etaeta,
       l_phi = l_phi, l_phiphi = l_phiphi, l_etaphi = l_etaphi)
}

#' Analytic gradient of the ZINB log-likelihood
#'
#' @inheritParams zinb_loglik_par
#' @return The gradient with respect to (beta, gamma, log(theta)).
#' @noRd
zinb_score <- function(par, X, Z, y, offset) {
  d <- zinb_deriv_pieces(par, X, Z, y, offset)
  c(crossprod(X, (1 - d$tau) * d$l_eta),
    crossprod(Z, d$tau - d$pi),
    sum((1 - d$tau) * d$l_phi))
}

#' Analytic Hessian of the ZINB log-likelihood
#'
#' @inheritParams zinb_loglik_par
#' @return The Hessian matrix with respect to (beta, gamma, log(theta)).
#' @noRd
zinb_hessian <- function(par, X, Z, y, offset) {
  d <- zinb_deriv_pieces(par, X, Z, y, offset)
  tt <- d$tau * (1 - d$tau)
  mt <- 1 - d$tau

  w_bb <- mt * d$l_etaeta + tt * d$l_eta^2
  w_gg <- tt - d$pi * (1 - d$pi)
  w_bg <- -tt * d$l_eta
  w_bp <- mt * d$l_etaphi + tt * d$l_eta * d$l_phi
  w_gp <- -tt * d$l_phi
  h_pp <- sum(mt * d$l_phiphi + tt * d$l_phi^2)

  H_bb <- crossprod(X, w_bb * X)
  H_gg <- crossprod(Z, w_gg * Z)
  H_bg <- crossprod(X, w_bg * Z)
  H_bp <- crossprod(X, w_bp)
  H_gp <- crossprod(Z, w_gp)

  rbind(cbind(H_bb,    H_bg,    H_bp),
        cbind(t(H_bg), H_gg,    H_gp),
        cbind(t(H_bp), t(H_gp), h_pp))
}
