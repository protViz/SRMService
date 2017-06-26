#' id lables for piw
#'
#' @export
#'
getIDLabels <- function(){
  .idlables <- c("Peptide.Sequence", "Protein.Name", "Precursor.Charge", "Product.Charge", "Fragment.Ion", "Isotope.Label")
}
#' piwots light
#' @param data data.frame containing columns : Peptide.Sequence ...
#'@export
piwotPiw <- function(data){
  lightPiw <- reshape2::dcast(data,
                              Peptide.Sequence + Protein.Name + Precursor.Charge +
                                Product.Charge + Fragment.Ion + Isotope.Label
                              ~ Replicate.Name,
                              value.var = "Area")
  idlables <- c("Peptide.Sequence", "Protein.Name", "Precursor.Charge", "Product.Charge", "Fragment.Ion", "Isotope.Label")
  return(lightPiw)
}

#' extracts Intensity columns
#' @param data data.frame returned by piwotPiw
#' @export
getIntensities <- function(data ){
  lightInt <- data[,7:ncol(data)]
  dum <- apply(data[ ,1:6],1, paste,collapse="_")
  rownames(lightInt) <- dum
  return(lightInt)
}


#' Make nice plots of transitions or peptides with correlations
#' @param dataX data.frame
#' @param main some name to plot
#' @param log log transform y axes
#' @param ylab label for y axes
#' @export
#'
plotNicely <- function(dataX, main="", log="", ylab="log(intensity)"){

  mat <- matrix(c(1,1,1,0,2,3), byrow=T, ncol=3)
  graphics::layout(mat, widths=c(5,2,2), heights=c(2,1))
  dataXt <- t(dataX)
  graphics::matplot(dataXt,type="b", main=main,lwd=1,lty=1, ylab="log(intensity)",las=2, xaxt = "n", log=log)
  graphics::axis(1, at = 1:nrow(dataXt), labels = rownames(dataXt), cex.axis = 0.7, las=2)
  graphics::legend("bottomleft", legend=rownames(dataX),col=1:5,lty=1 )
  nrow(dataX)
  if(nrow(dataX)>1){
    ordt <- (dataX)[order(apply(dataX,1,mean)),]
    dd <- stats::cor(t(ordt),use="pairwise.complete.obs", method = "spearman")
    imageWithLabelsNoLayout(dd,col=getBlueWhiteRed(),zlim=c(-1,1))
    imageColorscale(dd,col=getBlueWhiteRed(), zlim=c(-1,1))
    invisible(dd)
  }

}

#' Compute correlation matrix
#' @param dataX data.frame with transition intensities per peptide
#' @export
transitionCorrelations <- function(dataX){
  if(nrow(dataX) > 1){
    ordt <- (dataX)[order(apply(dataX,1,mean)),]
    dd <- stats::cor(t(ordt),use="pairwise.complete.obs", method = "spearman")
    return(dd)
  }else{
    message("Could not compute correlation, nr rows : " , nrow(dataX) )
  }

}






