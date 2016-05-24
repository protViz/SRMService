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

