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
                                     qValue = character(),

                                     # measurement levels
                                     hierarchy = list(),
                                     .hierarchy = list(),

                                     #h_fragmentId = character(),
                                     #h_precursorId = character(),
                                     #h_peptideId = character(),
                                     #h_proteinId = character(),


                                     isotopeLabel = character(),
                                     initialize = function(fileName="tmp"){
                                       self$fileName = fileName
                                     }
                                   )
)

ac <- AnalysisTableAnnotation$new()

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
  atable$isotopeLabel = isotopeLabel
  anaparam <- AnalysisParameters$new()
  configuration <- AnalysisConfiguration$new(atable, anaparam)
}

setupDataFrame <- function(data, configuration){
  required <- unlist(R6extractValues(configuration$table))
  required <- required[!grepl("^\\.", required)]
  longF <- dplyr::select(data,  required)

  for(i in 1:length(configuration$table$hierarchy))
  {
    longF <- unite(longF, UQ(sym(names(configuration$table$hierarchy)[i])), configuration$table$hierarchy[[i]],remove = FALSE)
  }

  attributes(longF)$configuration <- configuration
  return(longF)
}


#' Plot peptide and fragments
plotPeptides <- function(data,
                         proteinName,
                         sample,
                         intensity,
                         peptide,
                         fragment,
                         factor,
                         isotopeLabel
){
  print(as.list(match.call()))
  if(length(isotopeLabel)){
    formula <- paste(paste( isotopeLabel, collapse="+"), "~", paste(factor , collapse = "+"))
  }else{
    formula <- sprintf("~%s",factor)
  }

  p <- ggplot(data, aes_string(x = sample, y = intensity, group=fragment,  color= peptide))
  p <- p +  geom_point() + geom_line()
  p <- p + facet_grid(as.formula(formula), scales = "free_x"   )
  p <- p + scale_y_continuous(trans='log10')
  p <- p + ggtitle(proteinName) + theme_classic()
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top")
  return(p)
}

#' extracts the relevant information from the configuration to make the plot.
configurationPlotPeptides <- function(res, proteinName, configuration){
  hnames <- names(configuration$table$hierarchy)
  res <- plotPeptides(res, proteinName = proteinName,
                      sample = configuration$table$fileName,
                      intensity = configuration$table$startIntensity,
                      peptide = rev(hnames)[2],
                      fragment = rev(hnames)[1],
                      factor = configuration$table$factors[1],
                      isotopeLabel = configuration$table$isotopeLabel
  )
  return(res)
}





summarizeProtPepPrecursorFragCounts <- function(x, configuration){
  table <- configuration$table

  res <- list()
  for(i in 1:length(table$hierarchy)){
    name <- names(table$hierarchy)[i]
    values <- table$hierarchy[[i]]
    print(values)
    res[[name]] <- x %>% select(values) %>% distinct() %>% nrow()
  }
  res
}


#' Light only version.
summarizeProteins <- function(x, configuration ){
  hierarchy <- configuration$table$hierarchy

  precursorSum <- x %>% select(hierarchy$h_proteinId, hierarchy$h_precursorId, hierarchy$h_fragmentId) %>% distinct() %>%
    group_by(!!!(syms(c(hierarchy$h_proteinId, hierarchy$h_peptideId, hierarchy$h_precursorId)))) %>%
    summarize(nrFragments = n())


  peptideSum <- precursorSum %>% group_by(!!!(syms(c(hierarchy$h_proteinId,hierarchy$h_peptideId)))) %>%
    summarize(nrPrecursors = n(),
              minNrFragments = min(nrFragments),
              maxNrFragments = max(nrFragments))


  proteinSum <- peptideSum %>% group_by(!!!(syms(hierarchy$h_proteinId)))  %>%
    summarize(nrpeptides = n(),
              minNrPrecursors = min(nrPrecursors),
              maxNrPrecursors = max(nrPrecursors),
              maxNRFragments = max(maxNrFragments),
              minNrFragments= min(minNrFragments))

  proteinPeptide <- proteinSum %>% tidyr::unite(Precursors ,minNrPrecursors , maxNrPrecursors, sep="-", remove=FALSE)
  proteinPeptide <- proteinPeptide %>% tidyr::unite(Fragments ,minNrFragments , maxNRFragments, sep="-", remove=FALSE)
  return(proteinPeptide)
}

summarizeProteinsCounts <- function(x, configuration){
  hierarchy <- configuration$table$hierarchy

  precursor <- x %>% select(hierarchy$h_proteinId, hierarchy$h_precursorId, hierarchy$h_fragmentId) %>% distinct()
  x3<-precursor %>% group_by(!!!(syms(hierarchy$h_proteinId))) %>% summarize( nrPeptides =  n_distinct(!!!syms(hierarchy$h_peptideId)),
                                                                        nrPrecursors =  n_distinct(!!!syms(hierarchy$h_precursorId)),
                                                                        nrFragments =  n_distinct(!!!syms(hierarchy$h_fragmentId)))
}



