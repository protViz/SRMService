# Direct intensity manipulation ----
.setLargeQValuesToNA <- function(data,
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
#' @examples
#' analysis <- SRMService::spectronautDIAData250_analysis
#' config <- SRMService::spectronautDIAData250_config$clone(deep=TRUE)
#' res <- setLarge_Q_ValuesToNA(analysis, config)
setLarge_Q_ValuesToNA <- function(data, config, intensityNewName="IntensitiesWithNA"){
  data <- .setLargeQValuesToNA(data,
                               QValueColumn = config$table$qValue,
                               intensityOld = config$table$getWorkIntensity(),
                               thresholdQValue = config$parameter$maxQValue_Threshold,
                               intensityNew = intensityNewName
  )
  config$table$setWorkIntensity(intensityNewName)
  message("Column added and new WorkIntensity Set: ", intensityNewName)
  return(data)
}
#' DEPRECATED use \link{setLarge_Q_ValuesToNA}
#' @export
setIntensitiesToNA <- function(data, config, newcolname="IntensitiesWithNA"){
  warning("DEPRECATED use setLarge_Q_ValuesToNA instead")
  data <-setLarge_Q_ValuesToNA(data, config, intensityNewName=newcolname)
  return(data)
}

#' sets intensities smaller than threshold to NA
#' @export
#' @examples
#' analysis <- SRMService::spectronautDIAData250_analysis
#' config <- SRMService::spectronautDIAData250_config$clone(deep=TRUE)
#'
#' config$table$getWorkIntensity()
#'
#' config2 <- config$clone(deep=TRUE)
#' res1 <- setSmallIntensitiesToNA(analysis, config, threshold=1, intensityNewName = "Int1" )
#' res1000 <- setSmallIntensitiesToNA(analysis, config2, threshold=1000, intensityNewName = "Int1000" )
#' sum(is.na(res1[[config$table$getWorkIntensity()]])) < sum(is.na(res1000[[config2$table$getWorkIntensity()]]))
setSmallIntensitiesToNA <- function(data, config, threshold = 1 , intensityNewName ="IntensitiesWithNA"){
  resData <- data %>% mutate_at(vars(!!intensityNewName := config$table$getWorkIntensity()) ,
                                   function(x){ifelse(x < threshold, NA, x)} )
  config$table$setWorkIntensity(intensityNewName)
  return(resData)
}

#' Transform intensity
#' @export
#' @examples
#' library(tidyverse)
#' config <- SRMService::spectronautDIAData250_config$clone(deep=TRUE)
#' analysis <- SRMService::spectronautDIAData250_analysis
#' x <- transformIntensities(analysis, config, transform = log2)
#' stopifnot("log2_FG.Quantity" %in% colnames(x))
#' config <- SRMService::spectronautDIAData250_config$clone(deep=TRUE)
#' analysis <- SRMService::spectronautDIAData250_analysis
#' x <- transformIntensities(analysis, config, transform = asinh)
#' stopifnot("asinh_FG.Quantity" %in% colnames(x))
transformIntensities <- function(data,
                                 config,
                                 transformation,
                                 intesityNewName = NULL){
  x <- as.list( match.call() )
  if(is.null(intesityNewName)){
    newcol <- paste(as.character(x$transformation), config$table$getWorkIntensity(), sep="_")
  }else{
    newcol <- intesityNewName
  }
  data <- data %>% mutate_at(config$table$getWorkIntensity(), .funs = funs(!!sym(newcol) := transformation(.)))
  config$table$setWorkIntensity(newcol)
  print(x$transformation)
  if(grepl("log",as.character(x$transformation))){
    config$parameter$workingIntensityTransform = "log"
  }
  return(data)
}


# Summarize Q Values ----
#' Compute QValue summaries for each precursor
#' @export
#' @param data data
#' @param config configuration
#'
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
  message(glue("Columns added {QValueMin}, {QValueNR}"))
  return(data)
}




# Extact intensities in Wide format ----
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
#' @examples
#' res <- toWideConfig(sample_analysis, skylineconfig)
#' res <- toWideConfig(sample_analysis, skylineconfig, as.matrix = TRUE)
#' res <- scale(res)
#'
toWideConfig <- function(data, config , as.matrix = FALSE){
  res <- toWide( data, c(config$table$hierarchyKeys()[1],config$table$hierarchyKeys(TRUE)[1]) ,
                 config$table$sampleName ,
                 value = config$table$getWorkIntensity() )
  if(as.matrix){
    resMat <- as.matrix(select(res,-(1:2)))
    head(resMat)
    names <- res %>% select(1:2) %>% unite(precursor_id, 1,2, sep="~") %>% pull()
    rownames(resMat) <- names
    res <- resMat
  }
  return(res)
}

#' make it long
#' @export
#' @examples
#' conf <- skylineconfig$clone(deep = TRUE)
#' res <- toWideConfig(sample_analysis, conf, as.matrix = TRUE)
#' res <- scale(res)
#'
#' xx <- gatherItBack(res,"srm_intensityScaled",conf)
#' xx <- gatherItBack(res,"srm_intensityScaled",conf,sample_analysis)
#' conf$table$getWorkIntensity() == "srm_intensityScaled"
gatherItBack <- function(x,value,config,data = NULL){
  x <- dplyr::bind_cols(
    tibble::tibble("row.names" := rownames(x)),
    tibble::as_tibble(x)
  )
  x <- gather(x,key= !!config$table$sampleName, value = !!value, 2:ncol(x))
  x <- tidyr::separate(x, "row.names",  c(config$table$hierarchyKeys()[1],config$table$hierarchyKeys(TRUE)[1]), sep="~")
  if(!is.null(data)){
    x <- inner_join(data, x)
    config$table$setWorkIntensity(value)
  }
  return(x)
}

# Functions working on Matrices go Here ----
#' robust scale warpper
#' @export
robust_scale <- function(data){
  return(quantable::robustscale(data)$data)
}


#' apply Function To matrix
#' @export
#' @examples
#'
#' conf <- skylineconfig$clone(deep = TRUE)
#' res <- applyToIntensityMatrix(sample_analysis, conf, normalization = base::scale)
#' stopifnot("Area_base..scale" %in% colnames(res))
#' stopifnot("Area_base..scale" == conf$table$getWorkIntensity())
#'
#' res <- applyToIntensityMatrix(res, conf, normalization = robust_scale)
applyToIntensityMatrix <- function(data, config, normalization){
  x <- as.list( match.call() )
  colname <- make.names(paste(config$table$getWorkIntensity(), deparse(x$normalization), sep="_"))
  mat <- toWideConfig(data, config, as.matrix = TRUE)
  mat <- normalization(mat)
  data <- gatherItBack(mat, colname, config, data)
  return(data)
}



# Decorrelation analysis ----
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


# Missing Value imputation ----

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

  gatherItback <- function(x,config){
    x <- dplyr::bind_cols(
      tibble::tibble(!!config$table$hierarchyKeys(TRUE)[1] := rownames(x)),
      tibble::as_tibble(x)
    )
    gather(x,key= !!config$table$sampleName, value = "srm_ImputedIntensity", 2:ncol(x))
  }

  nestedX <- nestedX %>% dplyr::mutate(imputed = map(spreadMatrix, simpleImpute)) %>%
    dplyr::mutate(imputed = map(imputed, gatherItback, config))

  unnest_res <- nestedX %>% select(protein_Id, imputed) %>% unnest()
  qvalFiltX <- inner_join(x, unnest_res,
                          by=c(config$table$hierarchyKeys()[1], config$table$hierarchyKeys(TRUE)[1], config$table$sampleName) )
  config$table$setWorkIntensity("srm_ImputedIntensity")
  return(qvalFiltX)
}


#' Compute nr of B per A
#' @export
#' @examples
#' config <- skylineconfig$clone(deep=TRUE)
#' data <- sample_analysis
#' hierarchy <- config$table$hierarchyKeys()
#' res <- nr_B_in_A(data, hierarchy[1], hierarchy[2])
#' res %>% select(hierarchy[1],  nr_peptide_Id_by_protein_Id) %>%
#' distinct() %>% pull() %>% table()
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
  message("Column addded : ", c_name)
  return(data)
}

# Summarize Intensities by Intensity or NAs ----
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
#' @examples
#'
#' config <- spectronautDIAData250_config$clone(deep=T)
#' res <- setLarge_Q_ValuesToNA(spectronautDIAData250_analysis, config)
#' res <- rankPrecursorsByIntensity(res,config)
#' X <-res %>% select(c(config$table$hierarchyKeys(), srm_meanInt, srm_meanIntRank)) %>% distinct()
#' X %>% arrange(!!!syms(c(config$table$hierarchyKeys()[1], "srm_meanIntRank"  )))
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

#' aggregates top N intensities
#'
#' run \link{rankPrecursorsByIntensity} first
#' @export
#' @examples
#' config <- spectronautDIAData250_config$clone(deep=T)
#' res <- setLarge_Q_ValuesToNA(spectronautDIAData250_analysis, config)
#' res <- rankPrecursorsByIntensity(res,config)
#' aggregateTopNIntensities(res, config, N=3)
aggregateTopNIntensities <- function(data,config, N = 3){
  newcol <- "srm_sumTopInt"
  topInt <- data %>%
    dplyr::filter_at( "srm_meanIntRank", any_vars(. <= N)) %>%
    dplyr::group_by(!!!syms(c( config$table$hierarchyKeys()[1], config$table$sampleName)))
  sumNA <- function(x){sum(x, na.rm=TRUE)}
  sumTopInt <- topInt %>%
    dplyr::summarize( !!newcol := sumNA(!!sym(config$table$getWorkIntensity()))  )
  message("Column added : ", newcol)
  return(sumTopInt)
}

# Summarise NAs on lowest hierarchy ----

#' add number of NA's at lowest hierarchy to data
#' @export
#' @examples
#'
#' config <- spectronautDIAData250_config$clone(deep=T)
#' res <- setLarge_Q_ValuesToNA(spectronautDIAData250_analysis, config)
#' res <- summariseNAs(res,config)
#' x <- res %>%
#'   select(config$table$hierarchyKeys()[1], config$table$hierarchyKeys(T)[1], "srm_NrNAs") %>%
#'   distinct() %>% summarize(sum(srm_NrNAs)) %>% pull()
#' stopifnot(sum(is.na(res[[config$table$getWorkIntensity()[1]]])) == x)
summariseNAs <- function(data,
                         config){
  NrNAs = "srm_NrNAs"
  nNAs <- function(x){sum(is.na(x))}
  precursorID <- c(config$table$hierarchyKeys()[1], config$table$hierarchyKeys(TRUE)[1])
  workIntensity <- config$table$getWorkIntensity()

  lFG <- data %>% dplyr::group_by_at(precursorID)
  naSummaries <- lFG %>% dplyr::summarise_at( workIntensity, funs(!!NrNAs := nNAs(.)))
  data <- dplyr::inner_join(data, naSummaries , by=precursorID)
  message("Column added : ",NrNAs)
  return(data)
}


#' Ranks precursors by NAs (adds new column .NARank)
#' @export
#' @examples
#' config <- spectronautDIAData250_config$clone(deep=T)
#' res <- setLarge_Q_ValuesToNA(spectronautDIAData250_analysis, config)
#' res <- rankPrecursorsByNAs(res,config)
#' x <- res %>%
#'   select(config$table$hierarchyKeys()[1], config$table$hierarchyKeys(T)[1], "srm_NrNAs") %>%
#'   distinct() %>% summarize(sum(srm_NrNAs)) %>% pull()
#' stopifnot(sum(is.na(res[[config$table$getWorkIntensity()[1]]])) == x)
rankPrecursorsByNAs <- function(data, config){
  summaryColumn <- "srm_NrNAs"
  rankColumn <- "srm_NrNARank"
  data <- rankProteinPrecursors(data, config,
                               column = config$table$getWorkIntensity(),
                               fun = function(x){sum(is.na(x))},
                               summaryColumn = summaryColumn,
                               rankColumn = rankColumn,
                               rankFunction = function(x){min_rank(x)}
  )
  message("Added Columns :", summaryColumn, " ",  rankColumn)
  return(data)
}


#'


