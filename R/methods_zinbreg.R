#' Methods for zinbreg objects
#'
#' Standard S3 methods for models fitted by [zinbreg()].
#'
#' @param x,object An object of class `"zinbreg"` (or `"summary.zinbreg"`
#'   for the `print` method of a summary).
#' @param model For `coef` and `vcov`: `"full"` (both parts, the default),
#'   `"count"` or `"zero"`.
#' @param newdata An optional data frame with the covariates at which to
#'   predict. If omitted, predictions are for the data used in the fit.
#' @param type For `predict`: `"response"` (the mean \eqn{(1 - \pi)\mu},
#'   the default), `"count"` (the count-model mean \eqn{\mu}) or `"zero"`
#'   (the zero-inflation probability \eqn{\pi}). For `residuals`:
#'   `"response"` (`y - (1 - pi) mu`, the default) or `"pearson"` (the
#'   response residual divided by the ZINB standard deviation).
#' @param digits Number of significant digits to print.
#' @param ... Further arguments (currently ignored).
#'
#' @return
#' * `print` returns its argument invisibly.
#' * `summary` returns an object of class `"summary.zinbreg"` with separate
#'   coefficient tables (estimates, standard errors, z values, p-values)
#'   for the count and zero parts; the count table includes `Log(theta)`.
#' * `coef` returns the coefficients, named `count_*` and `zero_*` when
#'   `model = "full"`, and `vcov` their variance-covariance matrix.
#' * `logLik` returns an object of class `"logLik"` whose `df` counts both
#'   sets of coefficients plus `theta`, so that `AIC()` and `BIC()` work.
#' * `fitted`, `residuals` and `predict` return numeric vectors.
#' * `nobs` returns the number of observations.
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' d <- data.frame(x = rnorm(n), z = rnorm(n))
#' mu <- exp(1 + 0.5 * d$x)
#' p <- plogis(-1 + 0.8 * d$z)
#' d$y <- ifelse(runif(n) < p, 0, rnbinom(n, size = 2, mu = mu))
#' fit <- zinbreg(y ~ x | z, data = d)
#'
#' fit
#' summary(fit)
#' coef(fit)
#' coef(fit, model = "zero")
#' AIC(fit)
#' BIC(fit)
#' nd <- data.frame(x = c(-1, 0, 1), z = c(1, 0, -1))
#' predict(fit, newdata = nd)
#' predict(fit, newdata = nd, type = "count")
#' predict(fit, newdata = nd, type = "zero")
#'
#' @name zinbreg-methods
NULL

#' @rdname zinbreg-methods
#' @export
print.zinbreg <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nCount model coefficients (negative binomial with log link):\n")
  print(format(x$coefficients$count, digits = digits), quote = FALSE)
  cat("\nZero-inflation model coefficients (binomial with logit link):\n")
  print(format(x$coefficients$zero, digits = digits), quote = FALSE)
  cat("\nTheta: ", format(x$theta, digits = digits), "\n", sep = "")
  cat("Log-likelihood: ", format(x$loglik, digits = digits + 3), " on ",
      length(unlist(x$coefficients)) + 1, " Df\n\n", sep = "")
  invisible(x)
}

#' @rdname zinbreg-methods
#' @export
summary.zinbreg <- function(object, ...) {
  se <- sqrt(diag(object$vcov))
  make_table <- function(est, se) {
    z <- est / se
    cbind(Estimate = est, `Std. Error` = se, `z value` = z,
          `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  }

  b <- object$coefficients$count
  g <- object$coefficients$zero
  count_tab <- make_table(c(b, `Log(theta)` = log(object$theta)),
                          se[c(paste0("count_", names(b)), "Log(theta)")])
  zero_tab <- make_table(g, se[paste0("zero_", names(g))])
  rownames(count_tab) <- c(names(b), "Log(theta)")
  rownames(zero_tab) <- names(g)

  df <- length(b) + length(g) + 1
  out <- list(call = object$call,
              coefficients = list(count = count_tab, zero = zero_tab),
              theta = object$theta,
              SE.theta = object$SE.theta,
              loglik = object$loglik,
              df = df,
              aic = -2 * object$loglik + 2 * df,
              n = object$n,
              iterations = object$iterations,
              converged = object$converged)
  class(out) <- "summary.zinbreg"
  out
}

#' @rdname zinbreg-methods
#' @export
print.summary.zinbreg <- function(x,
                                  digits = max(3, getOption("digits") - 3),
                                  ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nCount model coefficients (negative binomial with log link):\n")
  stats::printCoefmat(x$coefficients$count, digits = digits,
                      signif.legend = FALSE)
  cat("\nZero-inflation model coefficients (binomial with logit link):\n")
  stats::printCoefmat(x$coefficients$zero, digits = digits,
                      signif.legend = TRUE)
  cat("\nTheta: ", format(x$theta, digits = digits),
      "  (SE ", format(x$SE.theta, digits = digits), ")\n", sep = "")
  cat("Log-likelihood: ", format(x$loglik, digits = digits + 3),
      " on ", x$df, " Df\n", sep = "")
  cat("AIC: ", format(x$aic, digits = digits + 3),
      "   Number of observations: ", x$n, "\n", sep = "")
  cat("Number of EM iterations: ", x$iterations,
      if (!x$converged) "  (did NOT converge)", "\n\n", sep = "")
  invisible(x)
}

#' @rdname zinbreg-methods
#' @export
coef.zinbreg <- function(object, model = c("full", "count", "zero"), ...) {
  model <- match.arg(model)
  b <- object$coefficients$count
  g <- object$coefficients$zero
  switch(model,
         count = b,
         zero = g,
         full = c(stats::setNames(b, paste0("count_", names(b))),
                  stats::setNames(g, paste0("zero_", names(g)))))
}

#' @rdname zinbreg-methods
#' @export
vcov.zinbreg <- function(object, model = c("full", "count", "zero"), ...) {
  model <- match.arg(model)
  nms <- names(coef(object, model = "full"))
  keep <- switch(model,
                 full = nms,
                 count = grep("^count_", nms, value = TRUE),
                 zero = grep("^zero_", nms, value = TRUE))
  V <- object$vcov[keep, keep, drop = FALSE]
  if (model != "full") {
    dimnames(V) <- list(sub("^(count|zero)_", "", keep),
                        sub("^(count|zero)_", "", keep))
  }
  V
}

#' @rdname zinbreg-methods
#' @export
logLik.zinbreg <- function(object, ...) {
  structure(object$loglik,
            df = length(unlist(object$coefficients)) + 1,
            nobs = object$n,
            class = "logLik")
}

#' @rdname zinbreg-methods
#' @export
nobs.zinbreg <- function(object, ...) {
  object$n
}

#' @rdname zinbreg-methods
#' @export
fitted.zinbreg <- function(object, ...) {
  object$fitted.values
}

#' @rdname zinbreg-methods
#' @export
residuals.zinbreg <- function(object, type = c("response", "pearson"),
                              ...) {
  type <- match.arg(type)
  res <- object$y - object$fitted.values
  if (type == "pearson") {
    # Var(Y) = (1 - pi) mu (1 + mu (pi + 1 / theta)) for the ZINB model.
    mu <- object$mu
    p <- object$pi
    res <- res / sqrt((1 - p) * mu * (1 + mu * (p + 1 / object$theta)))
  }
  res
}

#' @rdname zinbreg-methods
#' @export
predict.zinbreg <- function(object, newdata,
                            type = c("response", "count", "zero"), ...) {
  type <- match.arg(type)

  if (missing(newdata) || is.null(newdata)) {
    mu <- object$mu
    p <- object$pi
  } else {
    mf_old <- object$model
    mtX <- stats::delete.response(object$terms$count)
    mtZ <- object$terms$zero

    # Design matrices for the new data, with the factor levels and
    # contrasts of the original fit.
    X_old <- stats::model.matrix(object$terms$count, mf_old)
    Z_old <- stats::model.matrix(mtZ, mf_old)
    mfX <- stats::model.frame(mtX, newdata,
                              xlev = stats::.getXlevels(mtX, mf_old))
    mfZ <- stats::model.frame(mtZ, newdata,
                              xlev = stats::.getXlevels(mtZ, mf_old))
    X <- stats::model.matrix(mtX, mfX,
                             contrasts.arg = attr(X_old, "contrasts"))
    Z <- stats::model.matrix(mtZ, mfZ,
                             contrasts.arg = attr(Z_old, "contrasts"))

    off <- stats::model.offset(mfX)
    if (is.null(off)) off <- rep(0, nrow(X))
    if (!is.null(object$call$offset)) {
      off <- off + eval(object$call$offset, newdata)
    }

    mu <- exp(drop(X %*% object$coefficients$count) + off)
    p <- stats::plogis(drop(Z %*% object$coefficients$zero))
  }

  switch(type,
         response = (1 - p) * mu,
         count = mu,
         zero = p)
}
