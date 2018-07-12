library(R6)

# AnalysisParameters ----
#' Analysis parameters
#' @export
AnalysisParameters <- R6::R6Class("AnalysisParameters",
                                  public = list(
                                    maxQValue_Threshold  = 0.05,
                                    nrOfSigQvalues_Threshold = 5,
                                    qValThreshold = 0.01,
                                    minNumberOfQValues = 3,
                                    workingIntensityTransform = ""
                                  )
)

# AnalysisTableAnnotation ----
#' Create Annotation
#' @export
AnalysisTableAnnotation <- R6Class("AnalysisTableAnnotation",
                                   public = list(

                                     fileName = NULL,
                                     factors = list(), # ordering is important - first is considered the main

                                     sampleName = "sampleName",
                                     # measurement levels
                                     hierarchy = list(),
                                     retentionTime = NULL,
                                     isotopeLabel = character(),
                                     # do you want to model charge sequence etc?


                                     qValue = character(), # rename to score
                                     workIntensity = NULL, # could be list with names and functions
                                     startIntensity = NULL, # think of simplifying (use only workIntensity)

                                     initialize = function(){
                                     },
                                     setWorkIntensity = function(colName){
                                       self$workIntensity <- c(self$workIntensity, colName)
                                     },
                                     getWorkIntensity = function(){
                                       return(tail(self$workIntensity, n=1))
                                     },
                                     popWorkIntensity=function(){
                                        res <- self$workIntensity[length(self$workIntensity)]
                                        self$workIntensity <- self$workIntensity[-length(self$workIntensity)]
                                        return(res)
                                     },
                                     idRequired = function(){
                                       "Id Columns which must be in the input data frame"
                                       idVars <- c(
                                         self$fileName,
                                         unlist(self$factors),
                                         unlist(self$hierarchy),
                                         self$isotopeLabel
                                         )
                                       return(idVars)
                                     },
                                     hierarchyKeys = function(rev = FALSE){
                                       if(rev){
                                         return(rev(names(self$hierarchy)))
                                       }else{
                                         return(names(self$hierarchy))
                                       }
                                     },
                                     factorKeys = function(){
                                       return(names(self$factors))
                                     },
                                     idVars = function(){
                                       "Id Columns which must be in the output data frame"
                                       idVars <- c(
                                         self$fileName,
                                         names(self$factors),
                                         names(self$hierarchy),
                                         self$isotopeLabel,
                                         self$sampleName)
                                       return(idVars)
                                     },
                                     valueVars = function(){
                                       "Columns containing values"
                                       c(self$startIntensity, self$getWorkIntensity(), self$qValue)
                                     }
                                   )
)
# AnalysisConfiguration ----
#' Analysis Configuration
#' @export
AnalysisConfiguration <- R6Class("AnalysisConfiguration",
                                 public = list(
                                   table = NULL,
                                   parameter = NULL,
                                   initialize = function(analysisTableAnnotation, analysisParameter){
                                     self$table <- analysisTableAnnotation
                                     self$parameter <- analysisParameter
                                   }
                                 )
)

# Functions - Configuration ----
#' Helper function to extract all value slots in an R6 object
#' @export
R6extractValues <- function(r6class){
  tmp <- sapply(r6class, class)
  slots <- tmp[! tmp %in% c("environment", "function")]
  res <- list()
  for(i in names(slots)){
    if("R6" %in% class(r6class[[i]])){
      res[[i]]  <- R6extractValues(r6class[[i]])
    }else{
      res[[i]] <- r6class[[i]]
    }
  }
  return(res)
}
#' Deprecated
#' @export
setupDataFrame <- function(data, configuration ,sep="~"){
  warning("DEPRECATED replace with setup_analysis")
  res <- setup_analysis(data, configuration, sep)
  return(res)
}

#' Extracts columns relevant for a configuration from a data frame
#' @export
#' @examples
#'
#' skylineconfig <- craeteSkylineConfiguration(isotopeLabel="Isotope.Label.Type", qValue="Detection.Q.Value")
#' skylineconfig$table$factors[["Time"]] = "Sampling.Time.Point"
#' data(skylinePRMSampleData)
#'
#' sample_analysis <- setup_analysis(skylinePRMSampleData, skylineconfig)
#'
setup_analysis <- function(data, configuration ,sep="~"){
  table <- configuration$table
  for(i in 1:length(table$hierarchy))
  {
    data <- unite(data, UQ(sym(table$hierarchyKeys()[i])), table$hierarchy[[i]],remove = FALSE)
  }
  data <- select(data , -one_of(dplyr::setdiff(unlist(table$hierarchy), table$hierarchyKeys() )))

  for(i in 1:length(table$factors))
  {
    data <- unite(data, UQ(sym(table$factorKeys()[i])), table$factors[[i]],remove = FALSE)
  }

  sampleName <- table$sampleName

  if(!sampleName  %in% names(data)){
    message("creating sampleName")

    data <- data %>%  unite( UQ(sym( sampleName)) , unique(unlist(table$factors)), remove = TRUE ) %>%
      select(sampleName, table$fileName) %>% distinct() %>%
      mutate_at(sampleName, function(x){ x<- make.unique( x, sep=sep )}) %>%
      inner_join(data, by=table$fileName)
  } else{
    warning(sampleName, " already exists")
  }

  data <- data %>% select(-one_of(dplyr::setdiff(unlist(table$factors), table$factorKeys())))

  # Make implicit NA's explicit
  data <- data %>% select(c(configuration$table$idVars(),configuration$table$valueVars()))

  data <- complete( data , nesting(!!!syms(c(table$hierarchyKeys(), table$isotopeLabel))),
                    nesting(!!!syms(c( table$fileName , table$sampleName, table$factorKeys() ))))

  return( data )
}

# Functions - Plotting ----
#' Plot peptide and fragments
linePlotHierarchy_default <- function(data,
                                      proteinName,
                                      sample,
                                      intensity,
                                      peptide,
                                      fragment,
                                      factor,
                                      isotopeLabel,
                                      separate = FALSE,
                                      log_y=FALSE
){
  if(length(isotopeLabel)){
    if(separate){
      formula <- paste(paste( isotopeLabel, collapse="+"), "~", paste(factor , collapse = "+"))
      p <- ggplot(data, aes_string(x = sample,
                                   y = intensity,
                                   group=fragment,
                                   color= peptide
      ))
    }else{
      formula <- sprintf("~%s",factor)
      data <- unite(data, "fragment_label", fragment, isotopeLabel, remove = FALSE)
      p <- ggplot(data, aes_string(x = sample,
                                   y = intensity,
                                   group="fragment_label",
                                   color= peptide
      ))
    }
    p <- p +  geom_point(aes_string(shape= isotopeLabel)) + geom_line(aes_string(linetype = isotopeLabel))
  }else{
    formula <- sprintf("~%s",factor)
    p <- ggplot(data, aes_string(x = sample, y = intensity, group=fragment,  color= peptide))
    p <- p +  geom_point() + geom_line()

  }

  #p <- ggplot(data, aes_string(x = sample, y = intensity, group=fragment,  color= peptide, linetype = isotopeLabel))
  p <- p + facet_grid(as.formula(formula), scales = "free_x"   )
  p <- p + ggtitle(proteinName) + theme_classic()
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top")
  if(log_y){
    p <- p + scale_y_continuous(trans='log10')
  }
  return(p)
}

#' extracts the relevant information from the configuration to make the plot.
#' @export
#' @examples
#' library(SRMService)
#' conf <- SRMService::skylineconfig$clone(deep=TRUE)
#' xnested <- SRMService::sample_analysis %>%
#'  group_by_at(conf$table$hierarchyKeys()[1]) %>% tidyr::nest()
#'
#' SRMService::linePlotHierarchy_configuration(xnested$data[[1]], xnested$protein_Id[[1]],conf )
linePlotHierarchy_configuration <- function(res, proteinName, configuration, separate=FALSE){
  rev_hnames <- rev(names(configuration$table$hierarchy))
  res <- linePlotHierarchy_default(res, proteinName = proteinName,
                                   sample = configuration$table$sampleName,
                                   intensity = configuration$table$getWorkIntensity(),
                                   peptide = rev_hnames[2],
                                   fragment = rev_hnames[1],
                                   factor = names(configuration$table$factors)[1],
                                   isotopeLabel = configuration$table$isotopeLabel,
                                   separate = separate,
                                   log_y = (configuration$parameter$workingIntensityTransform != "log")
  )
  return(res)
}

#' add quantline to plot
#' @export
linePlotHierarchy_QuantLine <- function(p, data, aes_y,  configuration){
  table <- configuration$table
  p + geom_line(data=data,
                aes_string(x = table$sampleName , y = aes_y, group=1),
                size=1.3,
                color="black",
                linetype="dashed") +
    geom_point(data=data,
               aes_string(x = table$sampleName , y = aes_y, group=1), color="black", shape=10)
}

# Functions - summary ----

#' Count distinct elements for each level of hierarchy
#'
#' @export
#' @examples
#' library(SRMService)
#' skylineconfig <- craeteSkylineConfiguration(isotopeLabel="Isotope.Label.Type", qValue="Detection.Q.Value")
#' skylineconfig$table$factors[["Time"]] = "Sampling.Time.Point"
#' data(skylinePRMSampleData)
#'
#' sample_analysis <- setup_analysis(skylinePRMSampleData, skylineconfig)
#' hierarchyCounts(sample_analysis, skylineconfig)
hierarchyCounts <- function(x, configuration){
  hierarchy <- names( configuration$table$hierarchy )
  res <- x %>% group_by_at(configuration$table$isotopeLabel) %>% summarise_at( hierarchy, n_distinct )
  return(res)
}

#' Light only version.
#' Summarize Protein counts
#' @export
summarizeProteins <- function( x, configuration ){
  rev_hierarchy <- configuration$table$hierarchyKeys(TRUE)

  precursorSum <- x %>% select(rev_hierarchy) %>% distinct() %>%
    group_by_at(rev_hierarchy[-1]) %>%
    dplyr::summarize(nrFragments = n())

  peptideSum <- precursorSum %>% group_by_at(rev_hierarchy[-(1:2)]) %>%
    dplyr::summarize(nrPrecursors = n(),
              minNrFragments = min(nrFragments),
              maxNrFragments = max(nrFragments))


  proteinSum <- peptideSum %>% group_by_at(rev_hierarchy[-(1:3)])  %>%
    dplyr::summarize(nrpeptides = n(),
              minNrPrecursors = min(nrPrecursors),
              maxNrPrecursors = max(nrPrecursors),
              minNrFragments= min(minNrFragments),
              maxNrFragments = max(maxNrFragments)
              )
  proteinPeptide <- proteinSum %>% tidyr::unite(Precursors ,minNrPrecursors , maxNrPrecursors, sep="-", remove=FALSE)
  proteinPeptide <- proteinPeptide %>% tidyr::unite(Fragments ,minNrFragments , maxNrFragments, sep="-", remove=FALSE)
  return(proteinPeptide)
}
#' Summarize peptide Counts
#' @export
#' @examples
#' library(SRMService)
#' skylineconfig <- craeteSkylineConfiguration(isotopeLabel="Isotope.Label.Type", qValue="Detection.Q.Value")
#' skylineconfig$table$factors[["Time"]] = "Sampling.Time.Point"
#' data(skylinePRMSampleData)
#'
#' sample_analysis <- setup_analysis(skylinePRMSampleData, skylineconfig)
#' summarizeHierarchy(sample_analysis, skylineconfig)
#' summarizeHierarchy(sample_analysis, skylineconfig, level =2 )
#' summarizeHierarchy(sample_analysis, skylineconfig, level =3 )
#' summarizeHierarchy(sample_analysis, skylineconfig, level =4 )
summarizeHierarchy <- function(x, configuration, level = 1)
{
  hierarchy <- names(configuration$table$hierarchy)
  print(length(hierarchy))
  if(length(hierarchy) <= level){
    warning("There is less hierarchy levels than : ", level)
    return(NULL)
  }

  hierarchy <- hierarchy[level:length(hierarchy)]
  precursor <- x %>% select(hierarchy) %>% distinct()
  x3 <- precursor %>% group_by_at(hierarchy[1]) %>%
    dplyr::summarize_at( hierarchy[-1], n_distinct)
  return(x3)
}

# Functions - Missigness ----
#' compute missing statistics
#' @export
getMissingStats <- function(x, configuration, nrfactors = 1){
  table <- configuration$table
  factors <- head(table$factorKeys(), nrfactors)
  missingPrec <- x %>% group_by_at(c(factors,
                                     table$hierarchyKeys()[1],
                                     tail(table$hierarchyKeys(),1),
                                     table$isotopeLabel
  ))

  missingPrec <- missingPrec %>%
    dplyr::summarize(nrReplicates = n(), nrNAs = sum(is.na(!!sym(table$getWorkIntensity()))) ,
              meanArea = mean(!!sym(table$getWorkIntensity()), na.rm=TRUE)) %>%
    arrange(desc(nrNAs))
  missingPrec
}

#' Histogram summarizing missigness
#' @export
#' @examples
#' missignessHistogram(sample_analysis,skylineconfig)
#' setNa <- function(x){ifelse(x < 100, NA, x)}
#' sample_analysis %>% mutate(Area = setNa(Area)) -> sample_analysis
#' missignessHistogram(sample_analysis,skylineconfig)
missignessHistogram <- function(x, configuration, showempty = TRUE, nrfactors = 1){
  table <- configuration$table
  missingPrec <- getMissingStats(x, configuration,nrfactors)

  missingPrec <- missingPrec %>% ungroup()%>% dplyr::mutate(nrNAs = as.factor(nrNAs))
  if(showempty){
    if(configuration$parameter$workingIntensityTransform != "log")
    {
      missingPrec <- missingPrec %>% dplyr::mutate(meanArea = ifelse(is.na(meanArea),1,meanArea))
    }else{
      missingPrec <- missingPrec %>% dplyr::mutate(meanArea = ifelse(is.na(meanArea),-20,meanArea))
    }

  }

  factors <- head(table$factorKeys(), nrfactors)

  formula <- paste(table$isotopeLabel, "~", paste(factors, collapse = "+"))
  message(formula)

  p <- ggplot(missingPrec, aes(x = meanArea, fill = nrNAs, colour = nrNAs)) +
    geom_histogram(alpha = 0.2,position = "identity") +
    facet_grid(as.formula(formula)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  if(configuration$parameter$workingIntensityTransform != "log")
  {
    p <- p + scale_x_log10()
  }
  p
}

#' cumulative sums of missing
#' @export
#' @examples
#' setNa <- function(x){ifelse(x < 100, NA, x)}
#' sample_analysis %>% dplyr::mutate(Area = setNa(Area)) -> sample_analysis
#' res <- missingPerConditionCumsum(sample_analysis,skylineconfig)
#' names(res)
#' print(res$figure)
missingPerConditionCumsum <- function(x,configuration,nrfactors = 1){
  table <- configuration$table
  missingPrec <- getMissingStats(x, configuration,nrfactors)
  factors <- head(table$factorKeys(), nrfactors)

  xx <-missingPrec %>% group_by_at(c(table$isotopeLabel, factors,"nrNAs","nrReplicates")) %>%
    dplyr::summarize(nrTransitions =n())

  xxcs <-xx %>% group_by_at( c(table$isotopeLabel,factors)) %>% arrange(nrNAs) %>%
    dplyr::mutate(cs = cumsum(nrTransitions))
  res <- xxcs  %>% select(-nrTransitions)

  formula <- paste(table$isotopeLabel, "~", paste(factors, collapse = "+"))
  message(formula)
  p <- ggplot(res, aes(x= nrNAs, y = cs)) + geom_bar(stat="identity") +
    facet_grid(as.formula(formula))

  res <- res %>% spread("nrNAs","cs")
  return(list(data =res, figure=p))
}

#' Summarize missing in condtion as barplot
#' @export
#' @examples
#' setNa <- function(x){ifelse(x < 100, NA, x)}
#' sample_analysis %>% dplyr::mutate(Area = setNa(Area)) -> sample_analysis
#' res <- missingPerCondition(sample_analysis,skylineconfig)
#' names(res)
#' print(res$figure)
missingPerCondition <- function(x, configuration, nrfactors = 1){
  table <- configuration$table
  missingPrec <- getMissingStats(x, configuration, nrfactors)
  factors <- head(table$factorKeys(), nrfactors)

  xx <-missingPrec %>% group_by_at(c(table$isotopeLabel,
                                     factors,"nrNAs","nrReplicates")) %>%
    dplyr::summarize(nrTransitions =n())

  formula <- paste(table$isotopeLabel, "~", paste(factors, collapse = "+"))
  message(formula)

  p <- ggplot(xx, aes(x= nrNAs, y = nrTransitions)) + geom_bar(stat="identity")+
    facet_grid(as.formula(formula))

  xx <- xx %>% spread("nrNAs","nrTransitions")
  return(list(data = xx ,figure = p))
}

# Functions - Handling isotopes ----

#' spreads isotope label heavy light into two columns
#' @export
#' @examples
#' setNa <- function(x){ifelse(x < 100, NA, x)}
#' data(sample_analysis)
#' data(skylineconfig)
#' sample_analysis %>% dplyr::mutate(Area = setNa(Area)) -> sample_analysis
#' x<-spreadValueVarsIsotopeLabel(sample_analysis,skylineconfig)
#' head(x)
#'
#' x<-spreadValueVarsIsotopeLabel(sample_analysis_HL,skylineconfig_HL)
#' head(x[,5:ncol(x)])
spreadValueVarsIsotopeLabel <- function(resData, configuration){
  table <- configuration$table
  idVars <- table$idVars()
  resData2 <- resData %>% select(c(table$idVars(), table$valueVars()) )
  resData2 <- resData2 %>% gather(variable, value, - idVars  )
  resData2 <- resData2 %>%  unite(temp, table$isotopeLabel, variable )
  HLData <- resData2 %>% spread(temp,value)
  invisible(HLData)
}

# Computing protein Intensity summaries

.ExtractMatrix <- function(x){
  idx <- sapply(x,is.numeric)
  xmat <- as.matrix(x[,idx])
  rownames(xmat) <- x %>% select(which(!idx==TRUE)) %>% unite(x, sep="~") %>% pull(x)
  xmat
}

#' Extract intensity column
#' @export
#' @examples
#' library(dplyr)
#' xnested <- sample_analysis %>%
#'  group_by_at(skylineconfig$table$hierarchyKeys()[1]) %>%
#'  tidyr::nest()
#' xx <- extractIntensities(xnested$data[[1]],skylineconfig)
#' stopifnot(dim(xx)==c(104,22))
extractIntensities <- function(x, configuration){
  table <- configuration$table
  x <- x %>%
    select( c( table$sampleName,
               table$hierarchyKeys(TRUE)[1],
               table$getWorkIntensity()) ) %>%
    spread(table$sampleName, table$getWorkIntensity()) %>% .ExtractMatrix()
  return(x)
}

#' compute tukeys median polish from peptide or precursor intensities
#' @export
medpolishPly <- function(x){
  X <- medpolish(x,na.rm=TRUE, trace.iter=FALSE, maxiter = 10);
  res <- tibble("sampleName" = names(X$col) , medpolish = X$col + X$overall)
  res
}

#' realign data
#' @export
reestablishCondition <- function(data,
                                 medpolishRes,
                                 configuration
){
  table <- configuration$table
  xx <- data %>%  select(c(table$sampleName,
                           table$factorKeys(), table$fileName)) %>% distinct()
  res <- inner_join(xx,medpolishRes, by=table$sampleName)
  res
}


