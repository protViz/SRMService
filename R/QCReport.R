proteinColumns <- c("ProteinName","TopProteinName","nrPeptides")

#' Perform 2 group analysis with visualization
#' @export QCProteinReport
#' @exportClass QCProteinReport
#' @include eb.fit.R
#' @field proteinIntensity data.frame where colnames are Raw.File names, row.names are protein ID's and cells are protein abundances.
#' @field proteinAnnotation information about the proteins, nr of peptides etc.
#' @field nrPeptides min number of peptides per protein
#' @field maxNA maximum number of NA's
#' @field projectName name of project
#' @field experimentName name of experiment
#'
QCProteinReport <- setRefClass("QCProteinReport",
                            fields = list( proteinIntensity = "data.frame",
                                           proteinAnnotation = "data.frame",
                                           nrPeptides = "numeric",
                                           maxNA = "numeric",
                                           projectName = "character",
                                           experimentName = "character"
                            )
                            , methods = list(
                              setProteins = function( protein ){
                                "used to verify proteingroups structure and set members"
                                protein <- as.data.frame(protein)
                                stopifnot(proteinColumns %in% colnames(protein))
                                stopifnot(grep("Intensity\\." , colnames(protein)) > 0)
                                stopifnot(sum(duplicated(protein$TopProteinName))==0)
                                rownames(protein) <- protein$TopProteinName
                                protein <- protein[ protein$nrPeptides >= .self$nrPeptides, ]

                                .self$proteinIntensity <- protein[, grep("Intensity\\.",colnames(protein))]
                                colnames(.self$proteinIntensity) <- gsub("Intensity\\.","",colnames(.self$proteinIntensity))

                                # Sorts them in agreement with annotation_.
                                .self$proteinIntensity[.self$proteinIntensity==0] <- NA

                                nas <-.self$getNrNAs()
                                .self$proteinIntensity <- .self$proteinIntensity[nas <= maxNA,]
                                .self$proteinAnnotation <- protein[nas<=maxNA,proteinColumns]
                              },
                              initialize = function(
                                projectName,
                                experimentName="LFQ experiment",
                                maxNA=3,
                                nrPeptides = 2

                              ){
                                .self$projectName <- projectName
                                .self$experimentName <- experimentName
                                .self$nrPeptides <- nrPeptides
                                .self$maxNA <- maxNA
                              },
                              setMQProteinGroups = function(MQProteinGroups){
                                "set MQ protein groups table"
                                pint <- MQProteinGroups[,grep("Intensity\\.",colnames(MQProteinGroups))]
                                proteinTable <- data.frame(ProteinName = MQProteinGroups$Majority.protein.IDs,
                                                           TopProteinName = sapply(strsplit(MQProteinGroups$Majority.protein.IDs, split=";"),
                                                                                   function(x){x[1]}),
                                                           nrPeptides = MQProteinGroups$Peptides, pint, stringsAsFactors = F)
                                setProteins(proteinTable)
                              },
                              getNrNAs = function(){
                                'return number of NAs per protein'
                                return(quantable::rowNAs(.self$proteinIntensity))
                              },
                              getNormalized = function(){
                                quantable::robustscale(log2(.self$proteinIntensity))
                              }
                            )
)




