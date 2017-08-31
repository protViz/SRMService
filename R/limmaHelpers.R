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
#' get topTables for all coefficients from a limma lmfit ebayes result an merge them into single tibble.
#' @param lmfitebayes returned by lmFit
#' @param var since we work with tibbles what should be the colnames for row.names
#' @export
#'
getPVals <- function(lmfitebayes, var = "ProteinID"){
  res <- list()
  for(i in 1:length(colnames(lmfitebayes$coefficients)))
  {
    name <- colnames(lmfitebayes$coefficients)[i]
    res[[name]] <- data.frame(Condition = name, topTable(lmfitebayes, coef=name, number=Inf))
  }
  res <- lapply(res,tibble::rownames_to_column,var="ProteinID")
  res <- rbind.fill(res)
}

#' version of limma::contrast.fit which is able to handle NA's in one of the conditions.
#' @param fit returned by lmFit
#' @param cont returned by limma::makeContrasts
#' @return tibble with pValues and fold changes
#' @export
#' @examples
#' \dontrun{
#' cont <- limma::makeContrasts(contrasts = contrasts, levels =  levels)
#' fit <- lmFit(intmat , designMatrix)
#' res <- contrasts.fit.NA(fit,cont)
#' }
contrasts.fit.NA <- function(fit, cont){

  resl <- NULL
  for(i in 1:length(colnames(fit))){
    x <- colnames(fit)[i]
    #print(x)
    idx <- !grepl(x, colnames(cont))
    #print(colnames(cont)[idx])
    lmfit.cont <- contrasts.fit(fit[,-i], cont[-i,idx])
    lmfitebayes <- eBayes(lmfit.cont)
    resl[[i]] <- getPVals(lmfitebayes)
  }
  lmfit.cont <- contrasts.fit(fit, cont)
  lmfitebayes <- eBayes(lmfit.cont)

  resl <- c(resl, list(getPVals(lmfitebayes)))
  #return(resl)

  tmp <- rbind.fill(resl)
  nrNA <- apply(tmp,1,function(x){sum(is.na(x))})
  tmp <- tmp[nrNA == 0,]
  tmp <- tmp[!duplicated(subset(tmp, select =c("ProteinID","Condition"))),]

  return(dplyr::select(tmp, -adj.P.Val, -B))
}



