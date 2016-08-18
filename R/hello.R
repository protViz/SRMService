#' id lables for piw
#'
#' @export
#'
getIDLabels <- function(){
  .idlables <- c("Peptide.Sequence", "Protein.Name", "Precursor.Charge", "Product.Charge", "Fragment.Ion", "Isotope.Label")
}
#' piwots light
#'
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
#'
#' @export
getIntensities <- function(lightPiw ){
  lightInt <- lightPiw[,7:ncol(lightPiw)]
  dum <- apply(lightPiw[ ,1:6],1, paste,collapse="_")
  rownames(lightInt) <- dum
  return(lightInt)
}

#' get required columns for analysis
#'
#'@export
getRequiredColumns <- function(){
  cols <- c("Replicate.Name",
            "Peptide.Sequence",
            "Protein.Name",
            "Precursor.Charge",
            "Product.Charge",
            "Fragment.Ion",
            "Isotope.Label",
            "Precursor.Mz",
            "Product.Mz",
            "annotation_QValue",
            "Retention.Time",
            "Area",
            "Background")
  return(cols)
}

#' Make nice plots of transitions or peptides with correlations
#' @param dataX data.frame
#' @param main some name to plot
#' @export
#'
plotNicely <- function(dataX, main="", log="", ylab="log(intensity)"){

  mat <-matrix(c(1,1,1,0,2,3), byrow=T, ncol=3)
  layout(mat, widths=c(5,2,2), heights=c(2,1))
  dataXt <- t(dataX)
  matplot(dataXt,type="l", main=main,lwd=1,lty=1, ylab="log(intensity)",las=2, xaxt = "n", log=log)
  axis(1, at = 1:nrow(dataXt), labels = rownames(dataXt), cex.axis = 0.7, las=2)
  legend("bottomleft", legend=rownames(dataX),col=1:5,lty=1 )
  nrow(dataX)
  if(nrow(dataX)>1){
    ordt <- (dataX)[order(apply(dataX,1,mean)),]
    dd <- cor(t(ordt),use="pairwise.complete.obs", method = "spearman")
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
    dd <- cor(t(ordt),use="pairwise.complete.obs", method = "spearman")
    return(dd)
  }else{
    message("Could not compute correlation, nr rows : " , nrow(dataX) )
  }

}






