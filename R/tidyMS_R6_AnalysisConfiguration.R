library(R6)

AnalysisParameters <- R6::R6Class("AnalysisParameters",
                                  public = list(
                                    maxQValue_Threshold  = 0.05,
                                    nrOfSigQvalues_Threshold = 5,
                                    qValThreshold = 0.01,
                                    minNumberOfQValues = 3,
                                    workingIntensityTransform = ""
                                  )
)

#anaparam <- AnalysisParameters$new()
#anaparam$maxQValue_Threshold <- "tmp"

AnalysisTableAnnotation <- R6Class("AnalysisTableAnnotation",
                                   public = list(
                                     fileName = NULL,
                                     factors = list(), # ordering is important - first is considered the main
                                     sampleName = "sampleName",

                                     startIntensity = NULL,
                                     workIntensity = NULL, # could be list with names and functions


                                     retentionTime = NULL,
                                     qValue = character(),

                                     # measurement levels
                                     hierarchy = list(),

                                     isotopeLabel = character(),
                                     initialize = function(fileName="tmp"){
                                       self$fileName = fileName
                                     },
                                     idVars = function(){
                                       idVars <- c(
                                         self$fileName,
                                         names(self$factors),
                                         unlist(self$factors),
                                         names(self$hierarchy),
                                         unlist(self$hierarchy),
                                         self$isotopeLabel,
                                         self$sampleName)
                                       return(idVars)
                                     },
                                     idVars2 = function(){
                                       idVars <- c(
                                         self$fileName,
                                         names(self$factors),
                                         names(self$hierarchy),
                                         self$isotopeLabel,
                                         self$sampleName)
                                       return(idVars)
                                     },
                                     valueVars = function(){
                                       c(self$startIntensity, self$workIntensity, self$qValue)
                                     }
                                   )
)


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



#' helper function to extract all value slots in an R6 object
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


craeteSkylineConfiguration <- function(isotopeLabel="Isotope.Label", qValue="annotation_QValue"){
  atable <- AnalysisTableAnnotation$new()
  atable$fileName = "Replicate.Name"

  # measurement levels.
  atable$hierarchy[["protein_Id"]] <- "Protein.Name"
  atable$hierarchy[["peptide_Id"]] <- "Peptide.Sequence"
  atable$hierarchy[["precursor_Id"]] <-  c("Peptide.Sequence","Precursor.Charge")
  atable$hierarchy[["fragment_Id"]] <- c("Peptide.Sequence","Precursor.Charge","Fragment.Ion", "Product.Charge")


  #
  atable$qValue = qValue
  atable$startIntensity = "Area"
  atable$workIntensity = "Area"
  atable$isotopeLabel = isotopeLabel
  anaparam <- AnalysisParameters$new()
  configuration <- AnalysisConfiguration$new(atable, anaparam)
}

#' Extracts columns relevant for a configuration from a data frame
#'@export
setupDataFrame <- function(data, configuration ,sep="~"){
  table <- configuration$table

  for(i in 1:length(table$hierarchy))
  {
    data <- unite(data, UQ(sym(names(table$hierarchy)[i])), table$hierarchy[[i]],remove = FALSE)
  }
  data <- select(data , -one_of(setdiff(unlist(table$hierarchy), names(table$hierarchy))))

  for(i in 1:length(table$factors))
  {
    data <- unite(data, UQ(sym(names(table$factors)[i])), table$factors[[i]],remove = FALSE)
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

  data <- data %>% select(-one_of(setdiff(unlist(table$factors), names(table$factors))))

  # Make implicit NA's explicit

  data <- data %>% select(c(configuration$table$idVars2(),configuration$table$valueVars()))
  data <- complete( data , nesting(!!!syms(c(names(table$hierarchy), table$isotopeLabel))),
                    nesting(!!!syms(c( table$fileName , table$sampleName, names(table$factors) ))))

  attributes(data)$configuration <- configuration
  return( data )
}


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
linePlotHierarchy_configuration <- function(res, proteinName, configuration, separate=FALSE){
  rev_hnames <- rev(names(configuration$table$hierarchy))
  res <- linePlotHierarchy_default(res, proteinName = proteinName,
                                   sample = configuration$table$sampleName,
                                   intensity = configuration$table$workIntensity,
                                   peptide = rev_hnames[2],
                                   fragment = rev_hnames[1],
                                   factor = names(configuration$table$factors)[1],
                                   isotopeLabel = configuration$table$isotopeLabel,
                                   separate = separate,
                                   log_y = (configuration$parameter$workingIntensityTransform != "log")
  )
  return(res)
}


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


summarizeProtPepPrecursorFragCounts <- function(x, configuration){
  hierarchy <- names( configuration$table$hierarchy )
  res <- x %>% group_by_at(configuration$table$isotopeLabel) %>% summarise_at( hierarchy, n_distinct )
}




#' Light only version.
summarizeProteins <- function(x, configuration ){
  rev_hierarchy <- rev(names(configuration$table$hierarchy))
  print(rev_hierarchy)

  precursorSum <- x %>% select(rev_hierarchy) %>% distinct() %>%
    group_by_at(rev_hierarchy[-1]) %>%
    summarize(nrFragments = n())

  peptideSum <- precursorSum %>% group_by_at(rev_hierarchy[-(1:2)]) %>%
    summarize(nrPrecursors = n(),
              minNrFragments = min(nrFragments),
              maxNrFragments = max(nrFragments))


  proteinSum <- peptideSum %>% group_by_at(rev_hierarchy[-(1:3)])  %>%
    summarize(nrpeptides = n(),
              minNrPrecursors = min(nrPrecursors),
              maxNrPrecursors = max(nrPrecursors),
              maxNRFragments = max(maxNrFragments),
              minNrFragments= min(minNrFragments))
  proteinPeptide <- proteinSum %>% tidyr::unite(Precursors ,minNrPrecursors , maxNrPrecursors, sep="-", remove=FALSE)
  proteinPeptide <- proteinPeptide %>% tidyr::unite(Fragments ,minNrFragments , maxNRFragments, sep="-", remove=FALSE)
  return(proteinPeptide)
}

summarizeProteinsCounts <- function(x, configuration)
{
  rev_hierarchy <- rev( names( configuration$table$hierarchy ))
  hierarchy <- names(configuration$table$hierarchy)
  precursor <- x %>% select(rev_hierarchy) %>% distinct()
  x3<-precursor %>% group_by_at(hierarchy[1]) %>% summarize_at( hierarchy[-1], n_distinct)
  return(x3)
}


getMissingStats <- function(x, configuration, nrfactors = 1){
  table <- configuration$table
  factors <- head(names(table$factors), nrfactors)
  missingPrec <- x %>% group_by_at(c(factors,
                                     names(table$hierarchy)[1],
                                     tail(names(table$hierarchy),1),
                                     table$isotopeLabel
  )) %>%
    summarize(nrReplicates = n(), nrNAs = sum(is.na(!!sym(table$workIntensity))) ,
              meanArea = mean(!!sym(table$workIntensity), na.rm=TRUE)) %>%
    arrange(desc(nrNAs))
  missingPrec
}


missignessHistogram <- function(x, configuration, showempty = TRUE, nrfactors = 1){
  table <- configuration$table
  missingPrec <- getMissingStats(x, configuration,nrfactors)

  missingPrec <- missingPrec %>% ungroup()%>% mutate(nrNAs = as.factor(nrNAs))
  if(showempty){
    if(configuration$parameter$workingIntensityTransform != "log")
    {
      missingPrec <- missingPrec %>% mutate(meanArea = ifelse(is.na(meanArea),1,meanArea))
    }else{
      missingPrec <- missingPrec %>% mutate(meanArea = ifelse(is.na(meanArea),-20,meanArea))
    }

  }

  factors <- head(names(table$factors), nrfactors)

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


missingPerConditionCumsum <- function(x,configuration,nrfactors = 1){
  table <- configuration$table
  missingPrec <- getMissingStats(x, configuration,nrfactors)
  factors <- head(names(table$factors), nrfactors)

  xx <-missingPrec %>% group_by_at(c(table$isotopeLabel, factors,"nrNAs","nrReplicates")) %>%
    summarize(nrTransitions =n())

  xxcs <-xx %>% group_by_at( c(table$isotopeLabel,factors)) %>% arrange(nrNAs) %>% mutate(cs = cumsum(nrTransitions))
  res <- xxcs  %>% select(-nrTransitions)

  formula <- paste(table$isotopeLabel, "~", paste(factors, collapse = "+"))
  message(formula)
  p <- ggplot(res, aes(x= nrNAs, y = cs)) + geom_bar(stat="identity") +
    facet_grid(as.formula(formula))

  res <- res %>% spread("nrNAs","cs")
  return(list(data =res, figure=p))
}


missingPerCondition <- function(x, configuration, nrfactors = 1){
  table <- configuration$table
  missingPrec <- getMissingStats(x, configuration, nrfactors)
  factors <- head(names(table$factors), nrfactors)

  xx <-missingPrec %>% group_by_at(c(table$isotopeLabel,
                                     factors,"nrNAs","nrReplicates")) %>%
    summarize(nrTransitions =n())

  formula <- paste(table$isotopeLabel, "~", paste(factors, collapse = "+"))
  message(formula)

  p <- ggplot(xx, aes(x= nrNAs, y = nrTransitions)) + geom_bar(stat="identity")+
    facet_grid(as.formula(formula))

  xx <- xx %>% spread("nrNAs","nrTransitions")
  return(list(data = xx ,figure = p))
}


spreadValueVarsIsotopeLabel <- function(resData, configuration){
  table <- configuration$table
  idVars <- table$idVars()
  resData2 <- resData %>% select(c(table$idVars(), table$valueVars()) )
  resData2 <- resData2 %>% gather(variable, value, - idVars  )
  resData2 <- resData2 %>%  unite(temp, table$isotopeLabel, variable )
  HLData <- resData2 %>% spread(temp,value)
  invisible(HLData)
}



ExtractMatrix <- function(x){
  idx <- sapply(x,is.numeric)
  xmat <- as.matrix(x[,idx])
  rownames(xmat) <- x %>% select(which(!idx==TRUE)) %>% unite(x, sep="~") %>% pull(x)
  xmat
}


extractIntensities <- function(x, configuration){
  table <- configuration$table
  x <- x %>%
    select( c( table$sampleName,
               rev(names(table$hierarchy))[1],
               table$workIntensity) ) %>%
    spread(table$sampleName, table$workIntensity) %>% ExtractMatrix()
  return(x)
}


medpolishPly <- function(x, params){
  X <- medpolish(x,na.rm=TRUE, trace.iter=FALSE, maxiter = 10);
  res <- tibble("sampleName" = names(X$col) , medpolish = X$col + X$overall)
  res
}


reestablishCondition <- function(data,
                                 medpolishRes,
                                 configuration
){
  table <- configuration$table
  xx <- data %>%  select(c(table$sampleName,
                           names(table$factors), table$fileName)) %>% distinct()

  res <- inner_join(xx,medpolishRes, by=table$sampleName)
  res
}


