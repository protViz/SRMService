#' correlation preprocessing
#'
#' @export
#' @examples
#' rm(list=ls())
#' config <- spectronautDIAData250_config$clone(deep=T)
#' config$parameter$min_nr_of_notNA  <- 20
#' data <- spectronautDIAData250_analysis
#' res <- workflow_correlation_preprocessing(data,config)
workflow_correlation_preprocessing <- function(data, config, minCorrelation = 0.7){
  stat_input <- hierarchyCounts(data, config)

  data_NA <- setLarge_Q_ValuesToNA(data, config)
  data_NA <- summariseQValues(data_NA, config)

  data_NA_QVal <- data_NA %>% filter_at( "srm_QValueMin" , all_vars(. < config$parameter$qValThreshold )   )

  stat_qval <- hierarchyCounts(data_NA_QVal, config)

  # remove transitions with large numbers of NA's
  data_NA_QVal <- rankPrecursorsByNAs(data_NA_QVal, config)
  hierarchyCounts(data_NA_QVal, config)
  data_NA_QVal <- data_NA_QVal %>% dplyr::filter(srm_NrNotNAs > config$parameter$min_nr_of_notNA)

  stat_min_nr_of_notNA <- hierarchyCounts(data_NA_QVal, config)

  # remove single hit wonders
  data_NA_QVal <- nr_B_in_A(data_NA_QVal,config$table$hierarchyKeys()[1], config$table$hierarchyKeys()[2])
  c_name <- make_name(config$table$hierarchyKeys()[1], config$table$hierarchyKeys()[2])
  data_NA_QVal <- dplyr::filter(data_NA_QVal, !!sym(c_name) >= config$parameter$min_peptides_protein )

  stat_min_peptides_protein  <- hierarchyCounts(data_NA_QVal, config)

  # filter decorrelated.
  data_NA_QVal <- transformIntensities(data_NA_QVal, config, log2)
  data_NA_QVal <- markDecorrelated(data_NA_QVal, config, minCorrelation = minCorrelation)
  keepCorrelated <- dplyr::filter(data_NA_QVal, srm_decorelated == FALSE)

  stat_correlated  <- hierarchyCounts(keepCorrelated, config)

  # rank precursors by intensity
  keepCorrelated <- rankPrecursorsByIntensity(keepCorrelated, config)
  qvalFiltImputed <- impute_correlationBased(keepCorrelated, config)
  mean_na <- function(x){mean(x, na.rm = TRUE)}
  proteinIntensities <- aggregateTopNIntensities(qvalFiltImputed,config,func = mean_na,N=3)

  # collect stats
  stats <- list(stat_input = stat_input,
                stat_qval = stat_qval,
                stat_min_nr_of_notNA = stat_min_nr_of_notNA,
                stat_min_peptides_protein = stat_min_peptides_protein,
                stat_correlated = stat_correlated
                )
  x <- bind_rows(stats)
  stats <- add_column(x, processing = names(stats),.before = 1)

  return(list(data = proteinIntensities, stats = stats))
}
