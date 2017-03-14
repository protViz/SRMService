#' fit limma
#' @param dat data matrix
#' @param design design matrix
#' @export
#' @examples
#' data <- matrix(rnorm(120),ncol=6)
#' x <- c(rep("a",3),rep("b",3))
#' x <- model.matrix(~x)
#' eb.fit(data,x)
#'
eb.fit <- function(dat, design){
  fit <- limma::lmFit(dat, design)
  fit.eb <- limma::eBayes(fit)

  effectSize <- fit.eb$coefficients[, 2]
  df.r <- fit.eb$df.residual

  n <- nrow(dat)

  df.0 <- rep(fit.eb$df.prior, n)
  s2.0 <- rep(fit.eb$s2.prior, n)

  s2 <- (fit.eb$sigma)^2
  s2.post <- fit.eb$s2.post

  t.ord <- fit.eb$coefficients[, 2]/fit.eb$sigma/fit.eb$stdev.unscaled[, 2]
  t.mod <- fit.eb$t[, 2]

  p.ord <- 2*stats::pt(-abs(t.ord), fit.eb$df.residual)
  p.mod <- fit.eb$p.value[, 2]

  q.ord <- qvalue::qvalue(p.ord)$q
  q.mod <- qvalue::qvalue(p.mod)$q

  results.eb <- data.frame(effectSize, t.ord, t.mod, p.ord, p.mod, q.ord, q.mod, df.r, df.0, s2.0, s2, s2.post)
  results.eb <- results.eb[order(results.eb$p.mod), ]
  return(results.eb)
}
