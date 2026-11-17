test_that("fitting code uses no prohibited functions", {
  banned <- c("optim", "optimize", "optimise", "nlm", "nlminb", "uniroot",
              "glm", "glm.fit", "lm", "lm.fit", "lm.wfit",
              "glm.nb", "theta.ml", "zeroinfl",
              "optimx", "maxLik", "mle", "mle2")
  ns <- asNamespace("nbzinbSiam")
  used <- unlist(lapply(ls(ns, all.names = TRUE), function(nm) {
    obj <- get(nm, envir = ns)
    if (is.function(obj)) all.names(body(obj)) else character(0)
  }))
  used <- as.character(used)
  expect_length(intersect(used, banned), 0)
})

