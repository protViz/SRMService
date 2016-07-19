annotationColumns <- c("Raw.file","Condition", "BioReplicate", "Run","IsotopeLabelType")
proteinColumns <- c("ProteinName","TopProteinName","nrPeptides")

#'
#' @export Grp2Analysis
#' @exportClass Grp2Analysis
#'
#' @examples
#' Grp2Analysis()
Grp2Analysis <- setRefClass("Grp2Analysis",
                                        fields = list( proteinIntensity = "data.frame",
                                                       annotation = "data.frame",
                                                       proteinAnnotation = "data.frame",
                                                       nrPeptides = "numeric")
                                        , methods = list(
                                          setProteins = function(protein){
                                            "used to verify proteingroups structure and set members"
                                            stopifnot(proteinColumns %in% colnames(protein))
                                            stopifnot(grep("Intensity\\." , colnames(protein))>0)
                                            stopifnot(sum(duplicated(protein$TopProteinName))==0)
                                            rownames(protein) <- protein$TopProteinName
                                            .self$proteinIntensity <- protein[, grep("Intensity\\.",colnames(protein))]
                                            .self$proteinAnnotation <- protein[,proteinColumns]


                                          },
                                          initialize = function(
                                            annotation,
                                            nrPeptides = 2
                                          ){
                                            stopifnot(annotationColumns %in% colnames(annotation))
                                            .self$annotation <- annotation[order(annotation$Condition),]
                                          },
                                          setMQProteinGroups = function(MQProteinGroups){
                                            "set MQ protein groups table"
                                            pint <- MQProteinGroups[,grep("Intensity\\.",colnames(MQProteinGroups))]
                                            proteinTable <- data.frame(ProteinName = MQProteinGroups$Majority.protein.IDs,
                                                                       TopProteinName = sapply(strsplit(MQProteinGroups$Majority.protein.IDs, split=";"),
                                                                                               function(x){x[1]}),
                                                                       nrPeptides = MQProteinGroups$Peptides, pint, stringsAsFactors = F)
                                            setProteins(proteinTable)
                                          }
                                        )
)
