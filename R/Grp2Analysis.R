annotationColumns <- c("Raw.file","Condition", "BioReplicate", "Run","IsotopeLabelType")
proteinColumns <- c("ProteinName","TopProteinName","nrPeptides")

#' Perform 2 group analysis with visualization
#' @export Grp2Analysis
#' @exportClass Grp2Analysis
#'
Grp2Analysis <- setRefClass("Grp2Analysis",
                            fields = list( proteinIntensity = "data.frame",
                                           annotation = "data.frame",
                                           proteinAnnotation = "data.frame",
                                           nrPeptides = "numeric",
                                           maxNA = "numeric",
                                           conditions = "character",
                                           projectName = "character",
                                           experimentName = "character")
                            , methods = list(
                              setProteins = function(protein){
                                "used to verify proteingroups structure and set members"
                                stopifnot(proteinColumns %in% colnames(protein))
                                stopifnot(grep("Intensity\\." , colnames(protein))>0)
                                stopifnot(sum(duplicated(protein$TopProteinName))==0)
                                rownames(protein) <- protein$TopProteinName

                                protein <- protein[protein$nrPeptides >= .self$nrPeptides,]

                                .self$proteinIntensity <- protein[, grep("Intensity\\.",colnames(protein))]
                                colnames(.self$proteinIntensity) <- gsub("Intensity\\.","",colnames(.self$proteinIntensity))
                                stopifnot(colnames(.self$proteinIntensity) %in% .self$annotation$Raw.file)

                                .self$proteinIntensity <- .self$proteinIntensity[,.self$annotation$Raw.file]
                                .self$proteinIntensity[.self$proteinIntensity==0] <-NA

                                nas <-.self$getNrNAs()
                                .self$proteinIntensity <- .self$proteinIntensity[nas<=maxNA,]
                                .self$proteinAnnotation <- protein[nas<=maxNA,proteinColumns]

                              },
                              initialize = function(
                                annotation,
                                projectName,
                                experimentName="First experiment",
                                maxNA=3,
                                nrPeptides = 2
                              ){
                                .self$projectName <- projectName
                                .self$experimentName <- experimentName
                                stopifnot(annotationColumns %in% colnames(annotation))
                                .self$annotation <- annotation[order(annotation$Condition),]
                                .self$conditions <- unique(.self$annotation$Condition)
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
                                return(apply(.self$proteinIntensity, 1, function(x){sum(is.na(x))}))
                              },
                              getConditionData = function(condition){
                                'get intensities as matrix for single condition'
                                stopifnot(condition %in% .self$conditions)
                                fileID <-subset(.self$annotation, Condition == condition)$Raw.file
                                .self$proteinIntensity[, fileID]
                              },
                              getNormalized = function(){
                                quantable::robustscale(log2(.self$proteinIntensity))
                              },
                              getNormalizedConditionData = function(condition){
                                normalized <- .self$getNormalized()$data

                                stopifnot(condition %in% .self$conditions)
                                fileID <-subset(.self$annotation, Condition == condition)$Raw.file
                                normalized[, fileID]
                              },
                              getDesignMatrix = function(){
                                design <- model.matrix(~factor(.self$annotation$Condition))
                                return(design)
                              },
                              getPValues = function(){
                                return(eb.fit(.self$getNormalized()$data , .self$getDesignMatrix()))
                              }
                            )
)
