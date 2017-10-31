getConfig <- function(x){ attr(x,"configuration") }
setConfig <- function(x, config){ attributes(x)$configuration <- config; return(x) }

setupDataFrame <- function(data, configuration){
  required <- unique(unlist(configuration$required))
  longF <- select(data,  required)
  longF <- longF %>% unite(".PrecursorId", configuration$required$PrecursorId, remove = FALSE, sep=".")
  longF <- longF %>% unite(".SampleLabel", configuration$required$Factors, remove = FALSE, sep="_")
  configuration$.SampleLabel <- ".SampleLabel"
  configuration$.PrecursorId <- ".PrecursorId"
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


createWorkIntensities <- function(longF){
  config <- getConfig(longF)
  res <- longF
  threshold <- config$parameters$maxQValue_Threshold
  thresholdF <- function(x,y, threshold = 0.05){ ifelse(x < threshold, y, NA)}
  Vsyms <- rlang::syms(c(config$required$QValue,config$required$startIntensity ))
  res <- longF %>%
    mutate(!! config$workIntensity := thresholdF(!!!Vsyms))
  res <- setConfig(res, config)
  return(res)
}

m_inner_join <- function(x,y){
  config <- getConfig(x)
  x <- dplyr::inner_join(x,y)
  return(setConfig(x,config))
}

summariseQValues <- function(data){
  config <- getConfig(data)
  required <- config$required
  param <- config$parameters
  nthbestQValue <-  function(x,minNumberOfQValues){sort(x)[minNumberOfQValues]}
  npass <-  function(x,thresh = 0.05){sum(x < thresh)}
  lFG <- data %>%
    group_by_at(required$PrecursorId)
  qValueSummaries <- lFG %>% summarise_at(  c( required$QValue ),
                                            .funs = funs(".QValueMin"=nthbestQValue(.,param$minNumberOfQValues ),
                                                         ".QValueNR" =npass(., 0.05)
                                            ))

  return(qValueSummaries)
}


summariseNAs <- function(data){
  config <- getConfig(data)
  nNAs <- function(x ){sum(is.na(x))}
  lFG <- data %>%
    group_by_at(config$required$PrecursorId)

  naSummaries <- lFG %>% summarise_at( c(config$workIntensity), funs(".NrNAs" = nNAs(.)))
  longQNASummaries <- inner_join(data,naSummaries )


  config$.NrNAs=".NrNAs"
  longQNASummaries <- setConfig(longQNASummaries, config)

  return(longQNASummaries)
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

.removeDecorrelated <- function(ff, ValueCol , corThreshold = 0.7, tr = log2, .PrecursorId = ".PrecursorId" ){
  fx <-tr(ff[,ValueCol:ncol(ff)])
  rownames(fx) <- ff[,.PrecursorId]
  res <-.transitionCorrelations(fx, method="pearson")
  ff[!ff[,.PrecursorId] %in% .findDecorrelated(res,threshold = corThreshold),]
}


removeDecorellated <- function(data, minCorrelation = 0.65){
  config <- getConfig(data)
  required <- config$required

  tmp <- data %>%
    select_at(c(config$.PrecursorId, required$ProteinId, config$.SampleLabel,  config$workIntensity  )) %>%
    spread( key= config$.SampleLabel , value =  config$workIntensity )
  numericColumn <- length( c(config$.PrecursorId, required$ProteinId) ) + 1

  listProt <- plyr::dlply(tmp, required$ProteinId)

  intensitiesCorrelated <- plyr::llply(listProt, .removeDecorrelated,numericColumn , minCorrelation )
  intensitiesD <- plyr::ldply(intensitiesCorrelated, .id = required$ProteinId)
  return(intensitiesD)
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


## Ranks precursors by NAs (adds new column .NARank)
rankPrecursorsByNAs <- function(data, topNA = NULL){

  configuration <- getConfig(data)
  natmp <-data %>% group_by(!!!syms(c(configuration$.PrecursorId, configuration$required$ProteinId))) %>%
    summarise(.NrNAs = min(!!sym(configuration$.NrNAs)))
  x3 <- natmp %>% arrange(!!sym( configuration$required$ProteinId)) %>% group_by(!!sym( configuration$required$ProteinId))
  x3 <- x3 %>% mutate(.NARank = min_rank((.NrNAs)))
  if(!is.null(topNA)){
    x3 <- x3 %>% filter(.NARank < topNA)
  }
  configuration$.NARank <- ".NARank"
  data2 <- inner_join(data, x3)
  data <- setConfig(data2, configuration)
  return(data)
}


# impute missing
imputeMissing <- function(data){
  config <- getConfig(data)
  required <- config$required

  tmp <- data %>%
    select_at(c(config$.PrecursorId, required$ProteinId, config$.SampleLabel,  config$workIntensity  )) %>%
    spread( key= config$.SampleLabel , value =  config$workIntensity )
  numericColumn <- length( c(config$.PrecursorId, required$ProteinId) ) + 1

  listProt <- plyr::dlply(tmp, required$ProteinId)

  intensitiesImputed <- plyr::llply(listProt, imputef ,numericColumn )
  intensitiesImputed <- plyr::ldply(intensitiesImputed, .id = required$ProteinId)

  intensitiesImputedX <- gather(intensitiesImputed, !!sym(config$.SampleLabel),".ImputedIntensity",
                                -!!sym(config$.PrecursorId), -!!sym(required$ProteinId))

  head(intensitiesImputedX)
  data2 <- inner_join(data,intensitiesImputed)
  config$.ImputedIntensity <- ".ImputedIntensity"
  data2<-setConfig(data2,config)
  return(data2)
}





