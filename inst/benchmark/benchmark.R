# Benchmark nbzinbSiam against the reference implementations.
#
# For each data set (pscl::bioChemists and AER::NMES1988) and each model
# (NB and ZINB), this script fits the model with this package and with the
# reference function (MASS::glm.nb or pscl::zeroinfl), puts the estimates
# side by side, reports the maximum absolute differences (Table 4 of the
# project description) and times both fits.
#
# Run from the package's top-level folder, after installing the package:
#   source("inst/benchmark/benchmark.R")
# The results are saved to inst/benchmark/benchmark_results.rds, which the
# report reads.

library(nbzinbSiam)
library(MASS)
library(pscl)

data("bioChemists", package = "pscl")
data("NMES1988", package = "AER")

# Reference functions are run with tight tolerances, as the project
# description recommends.
nb_ref_control <- glm.control(epsilon = 1e-10, maxit = 100)
zinb_ref_control <- zeroinfl.control(reltol = 1e-12)

# Median elapsed time (seconds) over several runs of an expression.
time_it <- function(expr, times = 5) {
  expr <- substitute(expr)
  env <- parent.frame()
  median(replicate(times, system.time(eval(expr, env))[["elapsed"]]))
}

# ---------------------------------------------------------------------------
# NB: nbreg() vs MASS::glm.nb()
# ---------------------------------------------------------------------------
compare_nb <- function(formula, data, label) {
  fit <- nbreg(formula, data = data)
  ref <- glm.nb(formula, data = data, control = nb_ref_control)

  est <- data.frame(
    Parameter = c(names(coef(ref)), "theta"),
    nbreg = c(coef(fit), fit$theta),
    glm.nb = c(coef(ref), ref$theta),
    SE_nbreg = c(sqrt(diag(vcov(fit))), fit$SE.theta),
    SE_glm.nb = c(sqrt(diag(vcov(ref))), ref$SE.theta),
    row.names = NULL, check.names = FALSE
  )

  agreement <- data.frame(
    Data = label, Model = "NB",
    Coefficients = max(abs(coef(fit) - coef(ref))),
    SEs = max(abs(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(ref))))),
    Theta = abs(fit$theta - ref$theta),
    LogLik = abs(as.numeric(logLik(fit)) - as.numeric(logLik(ref))),
    LogLik_ours = as.numeric(logLik(fit)),
    LogLik_ref = as.numeric(logLik(ref))
  )

  timing <- data.frame(
    Data = label, Model = "NB",
    Ours = time_it(nbreg(formula, data = data)),
    Reference = time_it(glm.nb(formula, data = data,
                               control = nb_ref_control)),
    Iterations = fit$iterations
  )

  list(estimates = est, agreement = agreement, timing = timing)
}

# ---------------------------------------------------------------------------
# ZINB: zinbreg() vs pscl::zeroinfl()
# ---------------------------------------------------------------------------
compare_zinb <- function(formula, data, label) {
  fit <- zinbreg(formula, data = data)
  ref <- zeroinfl(formula, data = data, dist = "negbin",
                  control = zinb_ref_control)

  se_fit <- sqrt(diag(vcov(fit)))
  se_ref <- sqrt(diag(vcov(ref)))
  est <- data.frame(
    Parameter = c(names(coef(ref)), "Log(theta)"),
    zinbreg = c(coef(fit), log(fit$theta)),
    zeroinfl = c(coef(ref), log(ref$theta)),
    SE_zinbreg = c(se_fit, fit$SE.logtheta),
    SE_zeroinfl = c(se_ref, ref$SE.logtheta),
    row.names = NULL, check.names = FALSE
  )

  agreement <- data.frame(
    Data = label, Model = "ZINB",
    Coefficients = max(abs(coef(fit) - coef(ref))),
    SEs = max(abs(c(se_fit, fit$SE.logtheta) - c(se_ref, ref$SE.logtheta))),
    Theta = abs(log(fit$theta) - log(ref$theta)),
    LogLik = abs(as.numeric(logLik(fit)) - as.numeric(logLik(ref))),
    LogLik_ours = as.numeric(logLik(fit)),
    LogLik_ref = as.numeric(logLik(ref))
  )

  timing <- data.frame(
    Data = label, Model = "ZINB",
    Ours = time_it(zinbreg(formula, data = data), times = 3),
    Reference = time_it(zeroinfl(formula, data = data, dist = "negbin",
                                 control = zinb_ref_control), times = 3),
    Iterations = fit$iterations
  )

  list(estimates = est, agreement = agreement, timing = timing)
}

# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------
bio_nb <- art ~ fem + mar + kid5 + phd + ment
bio_zinb <- art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment

nmes_nb <- visits ~ hospital + health + chronic + gender + school + insurance
nmes_zinb <- visits ~ hospital + health + chronic + gender + school + insurance |
  hospital + chronic + gender + school + insurance

set.seed(651)
results <- list(
  bio_nb = compare_nb(bio_nb, bioChemists, "bioChemists"),
  bio_zinb = compare_zinb(bio_zinb, bioChemists, "bioChemists"),
  nmes_nb = compare_nb(nmes_nb, NMES1988, "NMES1988"),
  nmes_zinb = compare_zinb(nmes_zinb, NMES1988, "NMES1988")
)

agreement <- do.call(rbind, lapply(results, `[[`, "agreement"))
timing <- do.call(rbind, lapply(results, `[[`, "timing"))
rownames(agreement) <- rownames(timing) <- NULL

# Required tolerances (Table 4 of the project description).
tolerance <- data.frame(Model = c("NB", "ZINB"),
                        Coefficients = c(1e-5, 1e-3), SEs = c(1e-5, 1e-3),
                        Theta = c(1e-4, 1e-3), LogLik = c(1e-4, 1e-3))
tol <- tolerance[match(agreement$Model, tolerance$Model), -1]
agreement$Pass <- apply(agreement[, c("Coefficients", "SEs", "Theta",
                                      "LogLik")] <= tol, 1, all)

benchmark <- list(estimates = lapply(results, `[[`, "estimates"),
                  agreement = agreement, timing = timing,
                  tolerance = tolerance,
                  session = sessionInfo())

out_file <- file.path("inst", "benchmark", "benchmark_results.rds")
saveRDS(benchmark, out_file)

print(agreement, digits = 3)
print(timing, digits = 3)
cat("\nSaved to", out_file, "\n")
