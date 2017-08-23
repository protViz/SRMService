#' runs limma and computes qvalues
#' @param intMatrix intensity matrix
#' @param designMatrix design matrix
#' @return list with results
#' @export
#' @importFrom limma lmFit
#' @importFrom qvalue qvalue
#'
multigroupEBaysQvalue <- function(intMatrix, designMatrix){
  fit2 <- lmFit(intMatrix, designMatrix)
  fit2.eb <- eBayes(fit2)

  coefsigma <- with(fit2.eb, sweep(coefficients,1, sigma, "/"))
  t.ord <- with(fit2.eb, coefsigma/ stdev.unscaled )

  p.ord <- 2 * stats::pt(-abs(t.ord), fit2.eb$df.residual)
  p.mod <- fit2.eb$p.value

  q.ord <- qvalue::qvalue(p.ord)$q
  q.mod <- qvalue::qvalue(p.mod)$q

  coef <- (fit2.eb$coefficients)

  allElem <- list(intMatrix =intMatrix, effects = coef,
                  q.mod = q.mod,
                  q.ord = q.ord,
                  p.mod = p.mod,
                  p.ord = p.ord)
  allElem <- lapply(allElem, quantable::matrix_to_tibble)
  return(allElem)
}
#' merges results produced by multigroupEBaysQvalue
#' @param multigrpBayes result of multigroupEBaysQvalue function
#' @return data frame
#' @export
#' @importFrom plyr join_all
mergeLimmEBayesResult <- function(multigrpBayes){
  colnames(multigrpBayes$effects)[2:ncol(multigrpBayes$effects)] <- paste("Effect", colnames(multigrpBayes$effects)[2:ncol(multigrpBayes$effects)] , sep=".")
  colnames(multigrpBayes$q.mod)[2:ncol(multigrpBayes$q.mod)]  <- paste("q.mod", colnames(multigrpBayes$q.mod)[2:ncol(multigrpBayes$q.mod)], sep=".")
  colnames(multigrpBayes$q.ord)[2:ncol(multigrpBayes$q.ord)] <- paste("q.ord", colnames(multigrpBayes$q.ord)[2:ncol(multigrpBayes$q.ord)],  sep=".")
  colnames(multigrpBayes$p.ord)[2:ncol(multigrpBayes$p.ord)] <- paste("p.ord", colnames(multigrpBayes$p.ord)[2:ncol(multigrpBayes$p.ord)],  sep=".")
  colnames(multigrpBayes$p.mod)[2:ncol(multigrpBayes$p.mod)] <- paste("p.mod", colnames(multigrpBayes$p.mod)[2:ncol(multigrpBayes$p.mod)],  sep=".")
  colnames(multigrpBayes$intMatrix)[2:ncol(multigrpBayes$intMatrix)] <- paste("Intensity", colnames(multigrpBayes$intMatrix)[2:ncol(multigrpBayes$intMatrix)], sep=".")
  allData <- join_all(multigrpBayes, by="row.names")
  return(allData)
}
#' get pvalues for all coefficients from a limma lmfit ebayes result
#' @param lmfit.cont.ebayes
#' @param var = since we work with tibbles what should be the colnames for row.names
#' @export
#' @examples
getPVals <- function(lmfit.cont.ebayes, var = "ProteinID"){
  res <- list()
  for(i in 1:length(colnames(lmfit.cont.ebayes$coefficients)))
  {
    name <- colnames(lmfit.cont.ebayes$coefficients)[i]
    res[[name]] <- data.frame(Condition = name, topTable(lmfit.cont.ebayes, coef=name, number=Inf))
  }
  res <- lapply(res,tibble::rownames_to_column,var="ProteinID")
  res <- rbind.fill(res)
}

