library(R6)

AnalysisParameters <- R6::R6Class("AnalysisParameters",
                                  public = list(
                                    maxQValue_Threshold  = 0.05,
                                    nrOfSigQvalues_Threshold = 5,
                                    qValThreshold = 0.01,
                                    minNumberOfQValues = 3
                                  )
)

#anaparam <- AnalysisParameters$new()
#anaparam$maxQValue_Threshold <- "tmp"

AnalysisTableAnnotation <- R6Class("AnalysisTableAnnotation",
                                   public = list(
                                     #.workIntensity = NULL,
                                     fileName = NULL,
                                     factors = character(), # ordering is important - first is considered the main
                                     startIntensity = NULL,
                                     retentionTime = NULL,
                                     fragmentId = character(),
                                     precursorId = character(),
                                     peptideId = character(),
                                     proteinId = character(),
                                     qValue = character(),
                                     isotopeLabel = character(),
                                     .SampleLabel = character(),
                                     .fragmentId = character(),
                                     .precursorId = character(),
                                     initialize = function(fileName="tmp"){
                                       self$fileName = fileName
                                     }
                                   )
)

#ac <- AnalysisTableAnnotation$new()

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

R6extractValues <- function(r6class){
  tmp <- sapply(r6class, class)
  slots <- tmp[! tmp %in% c("environment", "function")]
  res <- list()
  for(i in names(slots)){
    res[[i]] <- r6class[[i]]
  }
  return(res)
}


craeteSkylineConfiguration <- function(isotopeLabel="Isotope.Label", qValue="annotation_QValue"){
  atable <- AnalysisTableAnnotation$new()
  atable$fileName = "Replicate.Name"
  atable$proteinId = "Protein.Name"
  atable$peptideId = "Peptide.Sequence"
  atable$precursorId = c("Peptide.Sequence","Precursor.Charge")
  atable$fragmentId = c("Peptide.Sequence","Precursor.Charge","Fragment.Ion", "Product.Charge")
  atable$qValue = qValue
  atable$startIntensity = "Area"
  atable$isotopeLabel = isotopeLabel
  anaparam <- AnalysisParameters$new()
  configuration <- AnalysisConfiguration$new(atable, anaparam)
}


#' Plot peptide and fragments
plotPeptides <- function(data,
                         proteinName,
                         x,
                         y,
                         group,
                         color,
                         factor,
                         isotopeLabel
){
  if(length(isotopeLabel)){
    formula <- paste(paste( isotopeLabel, collapse="+"), "~", paste(factor , collapse = "+"))
  }else{
    formula <- paste("~",factor, sep = "")
  }

  p <- ggplot(data, aes_string(x = x, y = y, group=".fragmentId", color="Peptide.Sequence" ))
  p <- p +  geom_point() + geom_line()
  p <- p + facet_grid(as.formula(formula), scales = "free_x"   )
  p <- p + scale_y_continuous(trans='log10')
  p <- p + ggtitle(proteinName) + theme_classic()
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top")
  return(p)
}

#' extracts the relevant information from the configuration to make the plot.
configurationPlotPeptides <- function(res, proteinName, configuration){
  res <- plotPeptides(res, proteinName = proteinName,
                      x = configuration$table$fileName,
                      y = configuration$table$startIntensity,
                      group = configuration$table$.fragmentId,
                      color = configuration$table$peptideId,
                      factor = configuration$table$factors[1],
                      isotopeLabel = configuration$table$isotopeLabel
  )
  return(res)
}


setupDataFrame <- function(data, configuration){
  required <- unlist(R6extractValues(configuration$table))
  required <- required[!grepl("^\\.", required)]
  print(required)
  longF <- dplyr::select(data,  required)

  longF <- longF %>% tidyr::unite(".fragmentId", configuration$table$fragmentId, remove = FALSE, sep=".")
  longF <- longF %>% tidyr::unite(".precursorId", configuration$table$precursorId, remove = FALSE, sep=".")
  longF <- longF %>% tidyr::unite(".SampleLabel", configuration$table$factors, remove = FALSE, sep="_")
  configuration$table$.SampleLabel = ".SampleLabel"
  configuration$table$.fragmentId = ".fragmentId"
  configuration$table$.precursorId = ".precursorId"
  attributes(longF)$configuration <- configuration
  return(longF)
}



summarizeProtPepPrecursorFragCounts <- function(x, configuration){
  table <- configuration$table
  fragment <- x %>% select(table$fragmentId) %>%
    distinct() %>% nrow()
  precursor <- x %>% select(table$precursorId) %>%
    distinct() %>% nrow()
  peptides <- x%>% select(table$peptideId) %>%
    distinct() %>% nrow()
  protein <- x %>% select(table$proteinId) %>%
    distinct() %>% nrow()
  c( nrproteins = protein , nrPeptides = peptides, nrPrecursor = precursor, nrFragments = fragment)
}


#' Light only version.
summarizeProteins <- function(x, configuration ){
  table <- configuration$table
  precursorSum <- x %>% select(table$proteinId, table$precursorId, table$fragmentId) %>% distinct() %>%
    group_by(!!!(syms(c(table$proteinId, table$peptideId, table$precursorId)))) %>%
    summarize(nrFragments = n())


  peptideSum <- precursorSum %>% group_by(!!!(syms(c(table$proteinId,table$peptideId)))) %>%
    summarize(nrPrecursors = n(),
              minNrFragments = min(nrFragments),
              maxNrFragments = max(nrFragments))


  proteinSum <- peptideSum %>% group_by(!!!(syms(table$proteinId)))  %>%
    summarize(nrpeptides = n(),
              minNrPrecursors = min(nrPrecursors),
              maxNrPrecursors = max(nrPrecursors),
              maxNRFragments = max(maxNrFragments),
              minNrFragments= min(minNrFragments))

  proteinPeptide <- proteinPeptide %>% tidyr::unite(Precursors ,minNrPrecursors , maxNrPrecursors, sep="-", remove=FALSE)
  proteinPeptide <- proteinPeptide %>% tidyr::unite(Fragments ,minNrFragments , maxNRFragments, sep="-", remove=FALSE)
  return(proteinPeptide)
}





