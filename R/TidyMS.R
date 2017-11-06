getConfig <- function(x){ attr(x,"configuration") }
setConfig <- function(x, config){ attributes(x)$configuration <- config; return(x) }
config_col_map <- function(x){ attr(x,"configuration")$required }
config_parameters <- function(x){ attr(x,"configuration")$parameters }

setupDataFrame <- function(data, configuration){
  required <- unique(unlist(configuration$required))
  longF <- select(data,  required)
  longF <- longF %>% unite(".PrecursorId", configuration$required$PrecursorId, remove = FALSE, sep=".")
  longF <- longF %>% unite(".SampleLabel", configuration$required$Factors, remove = FALSE, sep="_")
  configuration$.SampleLabel <- ".SampleLabel"
  configuration$.PrecursorId <- ".PrecursorId"

  message("Added Columns : .SampleLabel, .PrecursorId")
  attributes(longF)$configuration <- configuration
  return(longF)
}

summarizeCounts <- function(data){
  required <- getConfig(data)$required
  precursor <- unique(subset(data, select = required$PrecursorId))
  peptide <- unique(subset(data, select = required$PeptideId))
  protein <- unique(subset(data, select = required$ProteinId))
  list(nrprecursor = nrow(precursor), nrpeptide = nrow(peptide), nrproteins= nrow(protein ))
}


setIntensitiesToNA <- function(data,
                               threshold = config_parameters(data)$maxQValue_Threshold,
                               QValueColumn = config_col_map(data)$QValue,
                               intensityOld = config_col_map(data)$startIntensity,
                               intensityNew = getConfig(data)$workIntensity){
  config <- getConfig(data)

  thresholdF <- function(x,y, threshold = 0.05){ ifelse(x < threshold, y, NA)}
  data <- data %>%
    mutate(!! intensityNew := thresholdF(!!!syms(c(QValueColumn ,intensityOld )),threshold = threshold))


  data <- setConfig(data, config)
  return(data)
}

m_inner_join <- function(x,y){
  config <- getConfig(x)
  x <- dplyr::inner_join(x,y)
  return(setConfig(x,config))
}


#' Compute QValue summaries for each precursor
#' @param data 1
#' @param .PrecursorId 2
#' @param QValue 3
summariseQValues <- function(data,
                             precursorId = getConfig(data)$.PrecursorId,
                             QValue = config_col_map(data)$QValue,
                             maxQValThreshold = config_parameters(data)$maxQValue_Threshold,
                             minNumberOfQValues = config_parameters(data)$minNumberOfQValues
){
  config <- getConfig(data)

  .QValueMin <- ".QValueMin"
  .QValueNR <- ".QValueNR"

  nthbestQValue <-  function(x,minNumberOfQValues){sort(x)[minNumberOfQValues]}
  npass <-  function(x,thresh = maxQValThreshold){sum(x < thresh)}
  qValueSummaries <- data %>%
    group_by_at(precursorId) %>%
    summarise_at(  c( QValue ), .funs = funs(!!.QValueMin:=nthbestQValue(.,minNumberOfQValues ),
                                             !!.QValueNR :=npass(., maxQValThreshold)
    ))
  data <- inner_join(data, qValueSummaries)


  ## Add config...
  config$precursorStats[[.QValueMin]] <- .QValueMin
  config$precursorStats[[.QValueNR]] <- .QValueNR
  message("Added Columns:",.QValueMin,.QValueNR )
  data <- setConfig(data, config)

  return(data)
}


summariseNAs <- function(data,
                         precursorID = config_col_map(data)$PrecursorId,
                         workIntensity = getConfig(data)$workIntensity){
  config <- getConfig(data)

  .NrNAs = ".NrNAs"
  nNAs <- function(x){sum(is.na(x))}


  lFG <- data %>%
    group_by_at(precursorID)
  naSummaries <- lFG %>% summarise_at( workIntensity, funs(!!.NrNAs := nNAs(.)))
  data <- inner_join(data, naSummaries )


  ### add stuff to config.
  config$precursorStats[[.NrNAs]]=.NrNAs
  data <- setConfig(data, config)
  return(data)
}

### Filter using number of precursors
nrPrecursors <- function(data,
                         proteinID = getConfig(data)$required$ProteinId ,
                         precursorId = getConfig(data)$required$PrecursorId ){
  config <- getConfig(data)

  .nrPrecursors = ".nrPrecursors"

  tmp <- data %>%
    select_at(c(proteinID, precursorId)) %>%
    distinct() %>%
    group_by_at(proteinID) %>%
    summarise(!!.nrPrecursors:=n())
  data <- m_inner_join(data, tmp )

  ### add stuff to config.
  config$proteinStats[[.nrPrecursors]] <- .nrPrecursors
  data <- setConfig(data, config)
  return(data)
}


### Correlation Filtering
.transitionCorrelations <- function(dataX , method="spearman"){
  if(nrow(dataX) > 1){
    ordt <- (dataX)[order(apply(dataX,1,mean)),]
    dd <- stats::cor(t(ordt),use="pairwise.complete.obs", method = method)
    dd[is.na(dd)] <- -1
    return(dd)
  }else{
    message("Could not compute correlation, nr rows : " , nrow(dataX) )
  }
}

.findDecorrelated <- function(res, threshold = 0.65){
  nrtrans <- ncol(res)
  ids <- rowSums(res < threshold, na.rm = TRUE)
  names(which((nrtrans-1)== ids))
}

.removeDecorrelated <- function(ff, ValueCol , corThreshold = 0.7, tr = log2, .PrecursorId = ".PrecursorId" , .Decorrelated = ".Decorrelated"){
  fx <-tr(ff[,ValueCol:ncol(ff)])
  idcolumns <- ff[,1:(ValueCol-1)]
  rownames(fx) <- ff[,.PrecursorId]
  res <-.transitionCorrelations(fx, method="pearson")
  idcolumns[[.Decorrelated]] <- idcolumns[,.PrecursorId] %in% .findDecorrelated(res,threshold = corThreshold)
  return(idcolumns)
}

precursorIntensities2Wide <- function(data,
                                      precursorID = getConfig(data)$.PrecursorId,
                                      proteinID = config_col_map(data)$ProteinId,
                                      sampleLabel = getConfig(data)$.SampleLabel,
                                      value = config$workIntensity
){
  config <- getConfig(data)
  required <- config$required
  wide <- data %>%
    select_at(c(precursorID, proteinID, sampleLabel,  config$workIntensity  )) %>%
    spread( key= sampleLabel , value =  value )
  return(wide)
}


whichDecorellated <- function(data, minCorrelation = 0.65){
  config <- getConfig(data)
  required <- config$required
  newColumn <- ".Decorrelated"

  tmp <- precursorIntensities2Wide(data)
  numericColumn <- length( c(config$.PrecursorId, required$ProteinId) ) + 1

  listProt <- plyr::dlply(tmp, required$ProteinId)

  intensitiesCorrelated <- plyr::llply(listProt, .removeDecorrelated,numericColumn , minCorrelation )
  intensitiesD <- plyr::ldply(intensitiesCorrelated, .id = required$ProteinId)

  data <- inner_join(data, intensitiesD)
  config[[newColumn]]= newColumn
  data <- setConfig(data, config)
  message("Added Column: ", newColumn)
  return(data)
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


imputef <- function(xx,ValueCol){
  if(sum(is.na(xx)) == 0)
  {
    return(xx)
  }
  data <- data.frame(t(xx[,ValueCol:ncol(xx)]))
  # Why log transforming first
  data <- log2(data)

  c1<-simpleImpute(data)
  x <- data.frame(t(c1))

  x <- 2^x
  xx[,ValueCol:ncol(xx)] <- x
  if(sum(is.na(xx)) !=0){
    print("Still NA present")
  }
  return(xx)
}


# impute missing
imputeMissing <- function(data){
  imputedIntensityColname <- ".ImputedIntensity"
  config <- getConfig(data)
  required <- config$required

  tmp <- precursorIntensities2Wide(data)


  numericColumn <- length( c(config$.PrecursorId, required$ProteinId) ) + 1

  listProt <- plyr::dlply(tmp, required$ProteinId)

  intensitiesImputed <- plyr::llply(listProt, imputef ,numericColumn )
  intensitiesImputed <- plyr::ldply(intensitiesImputed, .id = required$ProteinId)

  intensitiesImputedX <- gather(intensitiesImputed, !!sym(config$.SampleLabel),!!sym(imputedIntensityColname),
                                -!!sym(config$.PrecursorId), -!!sym(required$ProteinId))

  head(intensitiesImputedX)
  data2 <- inner_join(data,intensitiesImputedX)
  config$.ImputedIntensity <- imputedIntensityColname
  message("Added Column:",imputedIntensityColname )
  data2<-setConfig(data2,config)
  return(data2)
}


rankProteinPrecursors <- function(data, top = NULL,
                                  column = configuration$workIntensity,
                                  fun = function(x){ mean(x, na.rm=TRUE)},
                                  summaryColumn = ".meanInt",
                                  rankColumn = ".meanIntRank",
                                  rankFunction = function(x){min_rank(desc(x))}
){
  configuration <- getConfig(data)
  summaryPerPrecursor <-data %>% group_by(!!!syms(c( configuration$.PrecursorId, configuration$required$ProteinId))) %>%
    summarise(!!summaryColumn := fun(!!sym(column)))

  groupedByProtein <- summaryPerPrecursor %>% arrange(!!sym( configuration$required$ProteinId)) %>% group_by(!!sym( configuration$required$ProteinId))
  rankedBySummary <- groupedByProtein %>% mutate(!!rankColumn := rankFunction(!!sym(summaryColumn)))

  data <- inner_join(data, rankedBySummary)
  return(data)
}


## Ranks precursors by NAs (adds new column .NARank)
rankPrecursorsPerProteinByNAs <- function(data){
  configuration <- getConfig(data)
  summaryColumn <- ".NrNAs"
  rankColumn <- ".NrNARank"

  data<- rankProteinPrecursors(data, column = configuration$precursorStats$.NrNAs,
                               fun = min,
                               summaryColumn = summaryColumn,
                               rankColumn = rankColumn,
                               rankFunction = function(x){min_rank(x)}
  )
  configuration[[summaryColumn]] <- summaryColumn
  configuration[[rankColumn]] <- rankColumn
  message("Added Columns :", summaryColumn,  rankColumn)
  data <- setConfig(data, configuration)
  return(data)
}


rankPrecursorsByIntensity <- function(data){
  configuration <- getConfig(data)
  summaryColumn <- ".meanInt"
  rankColumn <- ".meanIntRank"

  data<- rankProteinPrecursors(data, column = configuration$workIntensity,
                               fun = function(x){ mean(x, na.rm=TRUE)},
                               summaryColumn = summaryColumn,
                               rankColumn = rankColumn,
                               rankFunction = function(x){min_rank(desc(x))}
  )

  configuration[[summaryColumn]] <- summaryColumn
  configuration[[rankColumn]] <- rankColumn
  message("Added Columns :", summaryColumn, " ",  rankColumn)
  data <- setConfig(data, configuration)
  return(data)
}



aggregateTopNIntensities <- function(data, N = 3){
  config <- getConfig(data)
  required <- config$required

  topInt <- data %>% filter( .meanIntRank <= N) %>%
    group_by(!!!syms(c( required$ProteinId, config$.SampleLabel)))
  sumNA <- function(x){sum(x, na.rm=TRUE)}
  sumTopInt <-topInt %>% summarize( .sumTopInt =sumNA(!!sym(conf$workIntensity))  )

  config$proteinStats$.sumTopInt <- ".sumTopInt"
  data <- inner_join(data, sumTopInt)
  data <- setConfig(data, config)
  return(data)
}



