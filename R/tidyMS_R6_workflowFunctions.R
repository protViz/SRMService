#' correlation preprocessing
#'
#' @export
#' @examples
#' config <- spectronautDIAData250_config$clone(deep=T)
#'
#' data <- spectronautDIAData250_analysis
#' workflow_correlation_preprocessing(data,config)
workflow_correlation_preprocessing <- function(data, config){
  data_NA <- setLarge_Q_ValuesToNA(data, config)
  data_NA <- summariseQValues(data_NA, config)
  data_NA_QVal <- data_NA %>% filter_at( "srm_QValueMin" , all_vars(. < config$parameter$qValThreshold )   )

  # remove transitions with large numbers of NA's
  data_NA_QVal <- rankPrecursorsByNAs(data_NA_QVal, config)
  data_NA_QVal <- data_NA_QVal %>% dplyr::filter(srm_NrNotNAs > config$parameter$min_nr_of_notNA)

  # remove single hit wonders
  data_NA_QVal <- nr_B_in_A(data_NA_QVal,config$table$hierarchyKeys()[1], config$table$hierarchyKeys()[2])
  c_name <- make_name(config$table$hierarchyKeys()[1], config$table$hierarchyKeys()[2])
  data_NA_QVal <- dplyr::filter(data_NA_QVal, !!sym(c_name) >= config$parameter$min_peptides_protein )

  # filter decorrelated.
  data_NA_QVal <- transformIntensities(data_NA_QVal, config, log2)
  data_NA_QVal <- markDecorrelated(data_NA_QVal, config)
  removedCorrelated <- dplyr::filter(data_NA_QVal, srm_decorelated == FALSE)

  # rank precursors by intensity
  qvalFiltCorrSI <- rankPrecursorsByIntensity(qvalFiltCorr, config)
  qvalFiltImputed <- impute_correlationBased(qvalFiltCorrSI, config)
  proteinIntensities <- aggregateTopNIntensities(qvalFiltImputed,config,func = mean_na,N=3)
  return(proteinIntensities)
}
