#' List with peptide intensities
#'
#' @format A list of data frames
"correlatedPeptideList"

#' A data frame wich goes along with the \link{skylineconfig}.
#'
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{protein_Id}{protein id}
#'   \item{peptide_Id}{peptide id - stripped sequence}
#'   ...
#' }
#' @source \url{http://www.fgcz.ch/}
"sample_analysis"

#' A configuration which matches the \link{sample_analysis} data.
#'
#'
#' @format A AnalysisConfiguration R6 class
"skylineconfig"

#' Data frame which can be transformed into \link{sample_analysis}
#' by applying the \link{skylineconfig} using the function \link{setup_analysis}.
#'
#' @examples
#' data(skylineconfig)
#' data(skylinePRMSampleData)
#' x<-setup_analysis(skylinePRMSampleData,skylineconfig)
#' all.equal(x,sample_analysis)
#'
#' @format A AnalysisConfiguration R6 class
"skylinePRMSampleData"
