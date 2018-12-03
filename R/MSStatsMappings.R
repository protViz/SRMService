#' mapping of MSStats file names
#' @export
MSstatsMapping <- function(){
  MSstatsMapping <- list("ProteinName" = "Protein.Name",
                         "PeptideSequence" = "Peptide.Sequence",
                         "PrecursorCharge" = "Precursor.Charge",
                         "FragmentIon" = "Fragment.Ion",
                         "ProductCharge" = "Product.Charge",
                         "IsotopeLabelType"="Isotope.Label",
                         "Intensity"="Area",
                         "Condition"="Condition",
                         "BioReplicate"="Sample.Name",
                         "Run" = "Replicate.Name")
  return(MSstatsMapping)
}
