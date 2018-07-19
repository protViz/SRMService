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
                               QValueColumn = config$table$ident_qValue,
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
  message("Column added : ", newcol)
  if(grepl("log",as.character(x$transformation))){
    config$parameter$workingIntensityTransform = "log"
  }

  return(data)
}

#' visualize intensity distributions
#' @export
#' @examples
#' config <- skylineconfig$clone(deep=TRUE)
#' analysis <- transformIntensities(sample_analysis, config, log2)
#' plot_intensity_distribution_violin(analysis, config)
plot_intensity_distribution_violin <- function(data, config){
  ggplot(data, aes_string(x = config$table$sampleName, y = config$table$getWorkIntensity() )) +
    geom_violin() +
    theme(axis.text.x = element_text(angle=90))
}

#' visualize intensity distributions
#' @export
#' @examples
#' @rdname plot_intensity_distribution_violin
#' config <- skylineconfig$clone(deep=TRUE)
#' analysis <- transformIntensities(sample_analysis, config, log2)
#' plot_intensity_distribution_density(analysis, config)
plot_intensity_distribution_density <- function(data, config){
  ggplot(data, aes_string(x = config$table$getWorkIntensity(), colour = config$table$sampleName )) +
    geom_line(stat="density")
}

# Summarize Q Values ----
#' Compute QValue summaries for each precursor
#' @export
#' @param data data
#' @param config configuration
#' @examples
#' config <- skylineconfig$clone(deep=TRUE)
#' res <- summariseQValues(sample_analysis, config)
#' stopifnot(c("srm_QValueMin", "srm_QValueNR") %in% colnames(res))
summariseQValues <- function(data,
                             config
){
  QValueMin <- "srm_QValueMin"
  QValueNR <- "srm_QValueNR"

  precursorID <- config$table$hierarchyKeys(TRUE)[1]
  fileName <- config$table$fileName
  QValue  <- config$table$ident_qValue
  minNumberOfQValues <- config$parameter$minNumberOfQValues
  maxQValThreshold <- config$parameter$maxQValue_Threshold

  nthbestQValue <-  function(x,minNumberOfQValues){sort(x)[minNumberOfQValues]}
  npass <-  function(x,thresh = maxQValThreshold){sum(x < thresh)}

  qValueSummaries <- data %>%
    dplyr::select(fileName, precursorID, config$table$ident_qValue) %>%
    dplyr::group_by_at(precursorID) %>%
    dplyr::summarise_at(  c( QValue ), .funs = funs(!!QValueMin := nthbestQValue(.,minNumberOfQValues ),
                                                    !!QValueNR  := npass(., maxQValThreshold)
    ))
  print(colnames(qValueSummaries))
  data <- dplyr::inner_join(data, qValueSummaries, by=c(precursorID))
  message(glue::glue("Columns added {QValueMin}, {QValueNR}"))
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
markDecorrelated <- function(data , config, minCorrelation = 0.7){
  qvalFiltX <- data %>%  group_by_at(config$table$hierarchyKeys()[1]) %>% nest()
  qvalFiltX <- qvalFiltX %>%
    dplyr::mutate(spreadMatrix = map(data, extractIntensities, config))
  HLfigs2 <- qvalFiltX %>%
    dplyr::mutate(srmDecor = map(spreadMatrix, decorelatedPly, minCorrelation))
  unnest_res <- HLfigs2 %>%
    select(protein_Id, srmDecor) %>% unnest()
  qvalFiltX <- inner_join(data, unnest_res, by=c(config$table$hierarchyKeys()[1], config$table$hierarchyKeys(TRUE)[1]) )
  return(qvalFiltX)
}


# Missing Value imputation ----

simpleImpute <- function(data){
  m <-apply(data,2, mean, na.rm=TRUE )
  res <- sweep(data,2,m,"-")
  dim(data)
  dim(res)
  resMean <- apply(res, 1, mean, na.rm = TRUE)
  resid <- matrix(replicate(length(m),resMean), nrow=length(resMean))
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
  nestedX <- nestedX %>% dplyr::mutate(imputed = map(spreadMatrix, simpleImpute))
  nestedX <- nestedX %>% dplyr::mutate(imputed = map(imputed, gatherItback, config))

  unnest_res <- nestedX %>% select(protein_Id, imputed) %>% unnest()
  qvalFiltX <- inner_join(x, unnest_res,
                          by=c(config$table$hierarchyKeys()[1], config$table$hierarchyKeys(TRUE)[1], config$table$sampleName) )
  config$table$setWorkIntensity("srm_ImputedIntensity")
  return(qvalFiltX)
}

#' @export
make_name <- function(levelA, levelB, prefix="nr_"){
  c_name <- paste(prefix ,levelB,"_by_",levelA,sep="")
  return(c_name)
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
                      levelB, merge=TRUE){
  c_name <-make_name(levelA, levelB)
  tmp <- data %>%
    dplyr::select_at(c(levelA, levelB)) %>%
    dplyr::distinct() %>%
    dplyr::group_by_at(levelA) %>%
    dplyr::summarise(!!c_name:=n())
  if(!merge){
    return(tmp)
  }
  data <- dplyr::inner_join(data, tmp, by=levelA )
  message("Column added : ", c_name)
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

  message("Columns added:", summaryColumn, " ",  rankColumn)
  return(data)
}

#' aggregates top N intensities
#'
#' run \link{rankPrecursorsByIntensity} first
#' @export
#' @examples
#'
#' library(SRMService)
#' config <- spectronautDIAData250_config$clone(deep=T)
#' res <- setLarge_Q_ValuesToNA(spectronautDIAData250_analysis, config)
#' res <- rankPrecursorsByIntensity(res,config)
#' res %>% select(c(config$table$hierarchyKeys(),"srm_meanInt"  ,"srm_meanIntRank")) %>% distinct() %>% arrange(!!!syms(c(config$table$hierarchyKeys()[1],"srm_meanIntRank")))
#' mean_na <- function(x){mean(x, na.rm=TRUE)}
#' aggregateTopNIntensities(res, config, func = mean_na, N=3)
#'
aggregateTopNIntensities <- function(data , config, func, N){
  x <- as.list( match.call() )
  newcol <- make.names(glue::glue("srm_{deparse(x$func)}_{x$N}"))
  topInt <- data %>%
    dplyr::filter_at( "srm_meanIntRank", any_vars(. <= N)) %>%
    dplyr::group_by(!!!syms(c( config$table$hierarchyKeys()[1],
                               config$table$sampleName,
                               config$table$fileName,
                               config$table$factorKeys())))
  sumTopInt <- topInt %>%
    dplyr::summarize( !!newcol := func(!!sym(config$table$getWorkIntensity()))  )
  message("Column added : ", newcol)
  return(sumTopInt)
}

# Summarise NAs on lowest hierarchy ----

#' Ranks precursors by NAs (adds new column .NARank)
#' @export
#' @examples
#' config <- spectronautDIAData250_config$clone(deep=T)
#' res <- setLarge_Q_ValuesToNA(spectronautDIAData250_analysis, config)
#' res <- rankPrecursorsByNAs(res,config)
#' colnames(res)
#' x <- res %>%
#'   select(config$table$hierarchyKeys()[1], config$table$hierarchyKeys(T)[1], "srm_NrNotNAs") %>%
#'   distinct() %>% summarize(sum(srm_NrNotNAs)) %>% pull()
#' stopifnot(sum(!is.na(res[[config$table$getWorkIntensity()[1]]])) == x)
#' res %>% select(c(config$table$hierarchyKeys(),"srm_NrNotNAs"  ,"srm_NrNotNARank")) %>% distinct() %>% arrange(!!!syms(c(config$table$hierarchyKeys()[1],"srm_NrNotNARank")))
rankPrecursorsByNAs <- function(data, config){
  summaryColumn <- "srm_NrNotNAs"
  rankColumn <- "srm_NrNotNARank"
  data <- rankProteinPrecursors(data, config,
                                column = config$table$getWorkIntensity(),
                                fun = function(x){sum(!is.na(x))},
                                summaryColumn = summaryColumn,
                                rankColumn = rankColumn,
                                rankFunction = function(x){min_rank(desc(x))}
  )
  message("Columns added:", summaryColumn, " ",  rankColumn)
  return(data)
}
#'


