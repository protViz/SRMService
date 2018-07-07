# mark for deletion
### Correlation Filtering
.transitionCorrelations <- function(dataX , method="spearman"){
  if(nrow(dataX) > 1){
    ordt <- (dataX)[order(apply(dataX,1,mean)),]
    dd <- stats::cor(t(ordt),use="pairwise.complete.obs", method = method)
    dd[is.na(dd)] <- -1
    return(dd)
  }else{
    message("Could not compute correlation, nr rows : " , nrow(dataX) )
  }
}





.removeDecorrelated <- function(ff,
                                ValueCol ,
                                corThreshold = 0.7,
                                tr = log2,
                                .PrecursorId = ".PrecursorId" ,
                                .Decorrelated = ".Decorrelated"){
  fx <-tr(ff[,ValueCol:ncol(ff)])
  idcolumns <- ff[,1:(ValueCol-1)]
  rownames(fx) <- ff[,.PrecursorId]
  res <-.transitionCorrelations(fx, method="pearson")
  idcolumns[[.Decorrelated]] <- idcolumns[,.PrecursorId] %in% .findDecorrelated(res,threshold = corThreshold)
  return(idcolumns)
}
