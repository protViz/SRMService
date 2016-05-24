#' piwots light
#'
#'@export
piwotPiw <- function(light){
  lightPiw <- reshape2::dcast(light,Peptide.Sequence+Precursor.Charge+Protein.Name+ Fragment.Ion + Product.Charge ~ Replicate.Name, value.var = "Area")
  return(lightPiw)
}

#' extracts Intensity columns
#'
#' @export
getIntensities <- function(lightPiw ){
  lightInt <- lightPiw[,6:ncol(lightPiw)]
  dum <- apply(lightPiw[ ,c(1,2,3,4,5)],1, paste,collapse="_")
  rownames(lightInt) <- dum
  return(lightInt)
}


getRequiredColumns <- function(){
  cols <- c("Peptide.Sequence",
            "Protein.Name",
            "Replicate.Name",
            "Precursor.Charge",
            "Product.Charge",
            "Fragment.Ion",
            "annotation_QValue",
            "Isotope.Label",
            "Precursor.Mz",
            "Product.Mz",
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
plotNicely <- function(dataX, main=""){

  mat <-matrix(c(1,1,1,0,2,3), byrow=T, ncol=3)
  layout(mat, widths=c(5,2,2), heights=c(2,1))
  dataXt <- t(dataX)
  matplot(dataXt,type="l", main=main,lwd=1,lty=1, ylab="log(intensity)",las=2, xaxt = "n", log="y")

  axis(1, at = 1:nrow(dataXt), labels = rownames(dataXt), cex.axis = 0.7, las=2)
  legend("bottomleft", legend=rownames(dataX),col=1:5,lty=1 )
  nrow(dataX)
  if(nrow(dataX)>1){

    ordt <- (dataX)[order(apply(dataX,1,mean)),]

    dd <- cor(t(ordt),use="pairwise.complete.obs", method = "spearman")
    imageWithLabelsNoLayout(dd,col=getBlueWhiteRed(),zlim=c(-1,1))
    imageColorscale(dd,col=getBlueWhiteRed(), zlim=c(-1,1))
  }
}


