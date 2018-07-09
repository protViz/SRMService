.setIntensitiesToNA <- function(data,
                                QValueColumn,
                                intensityOld,
                                thresholdQValue = 0.05,
                                intensityNew = "IntensitiesWithNA"){

  thresholdF <- function(x,y, threshold = 0.05){ ifelse(x < threshold, y, NA)}
  data <- data %>%
    dplyr::mutate(!!intensityNew := thresholdF(!!!syms(c(QValueColumn ,intensityOld )),threshold = thresholdQValue))
  return(data)
}

#' sets intensities to NA if maxQValue_Threshold exceeded
#' @export
#' @return list with augmented data and updated config
setIntensitiesToNA <- function(data, config, newcolname="IntensitiesWithNA"){
  data <- .setIntensitiesToNA(data,
                              QValueColumn = config$table$qValue,
                              intensityOld = config$table$getWorkIntensity(),
                              thresholdQValue = config$parameter$maxQValue_Threshold,
                              intensityNew = newcolname
  )
  config$table$setWorkIntensity(newcolname)
  return(data)
}


m_inner_join <- function(x,y){
  config <- getConfig(x)
  x <- dplyr::inner_join(x,y)
  return(setConfig(x,config))
}


#' Compute QValue summaries for each precursor
#' @export
#' @param data data
#' @param config configuration
summariseQValues <- function(data,
                             config
){
  QValueMin <- "srm_QValueMin"
  QValueNR <- "srm_QValueNR"

  precursorID <- config$table$hierarchyKeys(TRUE)[1]
  fileName <- config$table$fileName
  QValue  <- config$table$qValue
  minNumberOfQValues <- config$parameter$minNumberOfQValues
  maxQValThreshold <- config$parameter$maxQValue_Threshold

  nthbestQValue <-  function(x,minNumberOfQValues){sort(x)[minNumberOfQValues]}
  npass <-  function(x,thresh = maxQValThreshold){sum(x < thresh)}

  qValueSummaries <- data %>%
    dplyr::select(fileName, precursorID, config$table$qValue) %>%
    dplyr::group_by_at(precursorID) %>%
    dplyr::summarise_at(  c( QValue ), .funs = funs(!!QValueMin := nthbestQValue(.,minNumberOfQValues ),
                                                    !!QValueNR  := npass(., maxQValThreshold)
    ))
  data <- dplyr::inner_join(data, qValueSummaries, by=c(precursorID))
  return(data)
}

#' add number of NA's per lowest hierarchy to data
#' @export
summariseNAs <- function(data,
                         config){
  NrNAs = "srm_NrNAs"
  nNAs <- function(x){sum(is.na(x))}
  precursorID <- config$table$hierarchyKeys(TRUE)[1]
  workIntensity <- config$table$getWorkIntensity()

  lFG <- data %>% dplyr::group_by_at(precursorID)
  naSummaries <- lFG %>% dplyr::summarise_at( workIntensity, funs(!!NrNAs := nNAs(.)))
  data <- dplyr::inner_join(data, naSummaries , by=precursorID)
  return(data)
}

#' Compute nr of B per A
#' @export
nr_B_in_A <- function(data,
                     levelA,
                     levelB){
  c_name <- paste("nr_",levelB,"_by_",levelA,sep="")
  tmp <- data %>%
    dplyr::select_at(c(levelA, levelB)) %>%
    dplyr::distinct() %>%
    dplyr::group_by_at(levelA) %>%
    dplyr::summarise(!!c_name:=n())
  data <- dplyr::inner_join(data, tmp, by=levelA )
  return(data)
}

#' transform long to wide
#' @export
toWide <- function(data,
                   rowIDs ,
                   columnLabels ,
                   value
){
  wide <- data %>%
    dplyr::select_at(c(rowIDs, columnLabels, value  ))
  wide <- wide %>%
    tidyr::spread( key= columnLabels , value =  value )
  return(wide)
}

#' transform long to wide
#' @export
toWideConfig <- function(data, config){
  return(toWide( data, config$table$hierarchyKeys(), config$table$sampleName , value = config$table$getWorkIntensity() ))
}

.findDecorrelated <- function(res, threshold = 0.65){
  if(is.null(res))
    return(NULL)
  nrtrans <- ncol(res)
  ids <- rowSums(res < threshold, na.rm = TRUE)
  names(which((nrtrans-1)== ids))
}

#' finds decorrelated measues
#' @export
decorelatedPly <- function(x, corThreshold = 0.7){
  res <- SRMService::transitionCorrelationsJack(x)
  decorelated <- .findDecorrelated(res,threshold = corThreshold)
  tibble(!!config$table$hierarchyKeys(TRUE)[1] := rownames(res), srm_decorelated = rownames(res) %in% decorelated)
}

#' marks decorrelated elements
#' @export
#' @section TODO: do something with warnings of type "the standard deviation is zero".
#' @section TODO: do investigate In max(x, na.rm = TRUE) : no non-missing arguments to max; returning -Inf
markDecorrelated <- function(x , config, minCorrelation = 0.7){
  x<-qvalFiltV
  qvalFiltX <- x %>%  group_by_at(config$table$hierarchyKeys()[1]) %>% nest()
  qvalFiltX <- qvalFiltX %>%
    dplyr::mutate(spreadMatrix = map(data, extractIntensities, config))
  HLfigs2 <- qvalFiltX %>%
    dplyr::mutate(srmDecor = map(spreadMatrix, decorelatedPly, minCorrelation))
  unnest_res <- HLfigs2 %>%
    select(protein_Id, srmDecor) %>% unnest()
  qvalFiltX <- inner_join(qvalFiltV, unnest_res, by=c(config$table$hierarchyKeys()[1], config$table$hierarchyKeys(TRUE)[1]) )
  return(qvalFiltX)
}


## Missing Value imputation

simpleImpute <- function(data){
  m <-apply(data,2, mean, na.rm=TRUE )
  res <- sweep(data,2,m,"-")
  resMean <- apply(res, 1, mean, na.rm = TRUE)
  resid <- replicate(length(m),resMean)
  imp <- sweep(resid,2,m,"+")
  res <- data
  res[is.na(res)] <- imp[is.na(res)]
  return(res)
}

#' imputation based on correlation assumption
#' @export
impute_correlationBased <- function(x , config){
  nestedX <- x %>%  group_by_at(config$table$hierarchyKeys()[1]) %>% nest()
  nestedX <- nestedX %>% dplyr::mutate(spreadMatrix = map(data, extractIntensities, config))

  spreadItback <- function(x,config){
    x <- dplyr::bind_cols(
      tibble::tibble(!!config$table$hierarchyKeys(TRUE)[1] := rownames(x)),
      tibble::as_tibble(x)
    )
    gather(x,key= !!config$table$sampleName, value = "srm_ImputedIntensity", 2:ncol(x))
  }

  nestedX <- nestedX %>% dplyr::mutate(imputed = map(spreadMatrix, simpleImpute)) %>%
    dplyr::mutate(imputed = map(imputed, spreadItback, config))

  unnest_res <- nestedX %>% select(protein_Id, imputed) %>% unnest()
  qvalFiltX <- inner_join(x, unnest_res,
                          by=c(config$table$hierarchyKeys()[1], config$table$hierarchyKeys(TRUE)[1], config$table$sampleName) )
  config$table$setWorkIntensity("srm_ImputedIntensity")
  return(qvalFiltX)
}


rankProteinPrecursors <- function(data,
                                  config,
                                  column = config$table$getWorkIntensity(),
                                  fun = function(x){ mean(x, na.rm=TRUE)},
                                  summaryColumn = "srm_meanInt",
                                  rankColumn = "srm_meanIntRank",
                                  rankFunction = function(x){min_rank(desc(x))}
){
  table <- config$table

  summaryPerPrecursor <-data %>%
    dplyr::group_by(!!!syms(table$hierarchyKeys())) %>%
    dplyr::summarise(!!summaryColumn := fun(!!sym(column)))

  groupedByProtein <- summaryPerPrecursor %>%
    dplyr::arrange(!!sym( table$hierarchyKeys()[1])) %>%
    dplyr::group_by(!!sym( table$hierarchyKeys()[1]))
  rankedBySummary <- groupedByProtein %>%
    dplyr::mutate(!!rankColumn := rankFunction(!!sym(summaryColumn)))

  data <- inner_join(data, rankedBySummary)
  return(data)
}

#' ranks precursor - peptide by intensity.
#'
#' @section TODO
#' @export
rankPrecursorsByIntensity <- function(data, config){
  summaryColumn <- "srm_meanInt"
  rankColumn <- "srm_meanIntRank"
  data<- rankProteinPrecursors(data, config, column = config$table$getWorkIntensity(),
                               fun = function(x){ mean(x, na.rm=TRUE)},
                               summaryColumn = summaryColumn,
                               rankColumn = rankColumn,
                               rankFunction = function(x){min_rank(desc(x))}
  )

  message("Added Columns :", summaryColumn, " ",  rankColumn)
  return(data)
}


#' Ranks precursors by NAs (adds new column .NARank)
#' @export
rankPrecursorsByNAs <- function(data, config){
  summaryColumn <- "srm_NrNAs"
  rankColumn <- "srm_NrNARank"
  data<- rankProteinPrecursors(data, config,
                               column = config$table$getWorkIntensity(),
                               fun = function(x){sum(is.na(x))},
                               summaryColumn = summaryColumn,
                               rankColumn = rankColumn,
                               rankFunction = function(x){min_rank(x)}
  )
  message("Added Columns :", summaryColumn, " ",  rankColumn)
  return(data)
}

#' aggregates top N intensities
#'
#' run \item{rankPrecursorsByIntensity} first
#' @export
aggregateTopNIntensities <- function(data,config, N = 3){
  topInt <- data %>%
    dplyr::filter_at( "srm_meanIntRank", any_vars(. <= N)) %>%
    dplyr::group_by(!!!syms(c( config$table$hierarchyKeys()[1], config$table$sampleName)))
  sumNA <- function(x){sum(x, na.rm=TRUE)}
  sumTopInt <- topInt %>%
    dplyr::summarize( srm_sumTopInt = sumNA(!!sym(config$table$getWorkIntensity()))  )
  #data <- inner_join(data, sumTopInt)
  return(sumTopInt)
}



