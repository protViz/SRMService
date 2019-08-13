#' Perform 2 group analysis with visualization
#' @export Grp2Analysis
#' @exportClass Grp2Analysis
#' @include eb.fit.R
#' @include RequiredColumns.R
#' @importFrom dplyr mutate
#' @field proteinIntensity data.frame where colnames are Raw.File names, row.names are protein ID's and cells are protein abundances.
#' @field annotation_ annotations data.frame with columns such as Raw.File, Condition, Run etc.
#' @field proteinAnnotation information about the proteins, nr of peptides etc.
#' @field nrPeptides min number of peptides per protein
#' @field maxNA maximum number of NA's
#' @field condition summary of conditions
#' @field projectName name of project
#' @field projectID name of experiment
#' @field pvalue pvalue threshold
#' @field qvalue qvalue threshold
#' @field pfoldchange foldchange threshold for p volcano plot
#' @field qfoldchange foldchange threshold for q volcano plot
#' @field reference document
#' @field removeDates strip the date from file name.
Grp2Analysis <- setRefClass("Grp2Analysis",
                            fields = list( proteinIntensity = "data.frame",
                                           annotation_ = "data.frame",
                                           proteinAnnotation = "data.frame",
                                           nrPeptides = "numeric",
                                           maxNA = "numeric",
                                           conditions = "character",
                                           projectName = "character",
                                           projectID = "character",
                                           workunitID = "character",
                                           pvalue= "numeric",
                                           qvalue= "numeric",
                                           pfoldchange = "numeric",
                                           qfoldchange = "numeric",
                                           reference = "character",
                                           removeDates= "logical",
                                           normalizationMethod = "character",
                                           housekeeper = "character",
                                           special = "character",
                                           modelMatrix = "matrix",
                                           numberOfProteinClusters = "numeric"
                            )
                            , methods = list(

                              initialize = function(
                                annotation,
                                projectName,
                                projectID="p1",
                                workunitID="wu1",
                                maxNA=3,
                                nrPeptides = 2,
                                reference = "Control",
                                annotationCol = annotationColumns,
                                removeDates = TRUE,
                                housekeeper = "",
                                normalizationMethod = "robustscale",
                                numberOfProteinClusters = 5
                              ){
                                if(!reference %in% unique(annotation$Condition)){
                                  stop("wrong reference :" , reference, " in conditions: ", paste(unique(annotation$Condition), collapse=","))
                                }
                                .self$projectName <- projectName
                                .self$projectID <- projectID
                                .self$workunitID <- workunitID
                                stopifnot(annotationCol %in% colnames(annotation))
                                annotation <- annotation[!is.na(annotation$Condition),]
                                .self$annotation_ <- annotation[order(annotation$Condition),]
                                .self$conditions <- as.character(unique(.self$annotation_$Condition))
                                .self$nrPeptides <- nrPeptides
                                .self$maxNA <- maxNA
                                .self$reference <- reference
                                .self$removeDates <- removeDates
                                .self$housekeeper <- housekeeper
                                .self$normalizationMethod <- normalizationMethod
                                .self$numberOfProteinClusters <- numberOfProteinClusters
                                setQValueThresholds()
                                setPValueThresholds()
                              },
                              setProteins = function(protein){
                                "used to verify proteingroups structure and set members"
                                protein <- as.data.frame(protein)
                                stopifnot(proteinColumns %in% colnames(protein))
                                stopifnot(grep("Intensity\\." , colnames(protein))>0)
                                stopifnot(sum(duplicated(protein$TopProteinName))==0)

                                rownames(protein) <- protein$TopProteinName
                                protein <- protein[protein$nrPeptides >= .self$nrPeptides,]

                                .self$proteinIntensity <- protein[, grep("Intensity\\.",colnames(protein))]
                                colnames(.self$proteinIntensity) <- gsub("Intensity\\.","",colnames(.self$proteinIntensity))

                                # Hack made for the FGCZ file conventions... Remove data from file start..

                                if(.self$removeDates == TRUE){
                                  colnames(.self$proteinIntensity) <- gsub("^[0-9]{8,8}_", "" ,colnames(.self$proteinIntensity))
                                  .self$annotation_$Raw.file <- gsub("^[0-9]{8,8}_", "" ,.self$annotation_$Raw.file)
                                }


                                stopifnot(.self$annotation_$Raw.file %in% colnames(.self$proteinIntensity))
                                .self$proteinIntensity <- .self$proteinIntensity[,.self$annotation_$Raw.file]
                                # Sorts them in agreement with annotation_.
                                .self$proteinIntensity[.self$proteinIntensity==0] <- NA

                                nas <-.self$getNrNAs()
                                .self$proteinIntensity <- .self$proteinIntensity[nas <= maxNA,]
                                .self$proteinAnnotation <- protein[nas<=maxNA,proteinColumns]

                              },
                              setQValueThresholds = function(qvalue = 0.05, qfoldchange=0.1){
                                .self$qvalue = qvalue
                                .self$qfoldchange = qfoldchange
                              },
                              setPValueThresholds = function(pvalue = 0.01, pfoldchange=0.5){
                                .self$pvalue= pvalue
                                .self$pfoldchange = pfoldchange
                              },
                              setMQProteinGroups = function(MQProteinGroups){
                                "set MQ protein groups table"
                                pint <- MQProteinGroups[,grep("Intensity\\.",colnames(MQProteinGroups))]
                                proteinTable <- data.frame(ProteinName = MQProteinGroups$Majority.protein.IDs,
                                                           TopProteinName = sapply(strsplit(MQProteinGroups$Majority.protein.IDs, split=";"),
                                                                                   function(x){x[1]}),
                                                           nrPeptides = MQProteinGroups$Peptides,
                                                           Fasta.headers = MQProteinGroups$Fasta.headers,
                                                           pint, stringsAsFactors = F)
                                setProteins(proteinTable)
                              },
                              getNrNAs = function(){
                                'return number of NAs per protein'
                                return(apply(.self$proteinIntensity, 1, function(x){sum(is.na(x))}))
                              },
                              getConditionData = function(condition){
                                'get intensities as matrix for single condition'
                                stopifnot(condition %in% .self$conditions)
                                fileID <- as.character(subset(.self$annotation_, Condition == condition)$Raw.file)
                                .self$proteinIntensity[, fileID]
                              },
                              setNormalizationMethod = function(normalizationMethod = "robustscale", housekeeper = ""){
                                'set the normalization parameters'
                                .self$housekeeper <- housekeeper
                                .self$normalizationMethod <- normalizationMethod
                              },
                              getNormalized = function(){
                                'return normalized data'
                                normMethod <- .self$normalizationMethod
                                if(normMethod == "robustscale"){
                                  quantable::robustscale(log2(.self$proteinIntensity))
                                }else if(normMethod == "vsn"){
                                  message("not implemented")
                                }else if(normMethod == "none"){
                                  dumm <-quantable::robustscale(log2(.self$proteinIntensity))
                                  list(data=log2(.self$proteinIntensity), medians=dumm$medians, mads=dumm$mads)
                                }else if(normMethod == "housekeeper"){
                                  message("not implemented")
                                }else{
                                  stop(normMethod, "is not a valid normalizaiton method")
                                }

                              },
                              getNormalizedVSN = function(){
                                "perform variance stabilizing normalizaiton."

                              },
                              getNormalizedConditionData = function(condition){
                                normalized <- .self$getNormalized()$data
                                stopifnot(condition %in% .self$conditions)
                                fileID <- as.character(subset(.self$annotation_, Condition == condition)$Raw.file)
                                normalized[, fileID]
                              },
                              resetModelMatrix = function(){
                                .self$modelMatrix <- matrix(ncol=0, nrow=0)
                              },
                              setModelMatrix = function(modelMatrix){
                                .self$modelMatrix <- modelMatrix
                              },
                              getDesignMatrix = function(){
                                if(isTRUE(all.equal(dim(.self$modelMatrix), c(0,0)))){
                                  # hack for putting reference into denomintor.
                                  design <- forcats::fct_relevel(.self$annotation_$Condition, .self$reference)
                                  .self$modelMatrix <- model.matrix(~design)
                                }
                                return(.self$modelMatrix)
                              },
                              #getPValues = function(){
                              #  res <- eb.fit(.self$getNormalized()$data , .self$getDesignMatrix())
                              #  res <- data.frame(proteinID = rownames(res), log2FC = res$effectSize, res)
                              #  return(res)
                              #},
                              getFit = function(){
                                fit <- limma::lmFit(.self$getNormalized()$data, .self$getDesignMatrix())
                                fit.eb <- limma::eBayes(fit)
                                return(list(fit= fit, fit.eb = fit.eb))
                              },
                              getModPValuesCI = function(){
                                "get Pvalues and confidence intervals"
                                fit <- limma::lmFit(.self$getNormalized()$data, .self$getDesignMatrix())
                                fit.eb <- limma::eBayes(fit)
                                res <- topTable(fit.eb, coef=2, number=Inf,confint = TRUE)
                                res <- data.frame(proteinID = rownames(res), log2FC = res$logFC, res)
                                res$logFC <- NULL
                                return(res)
                              },
                              getAnnotation = function(){
                                return(.self$annotation_)
                              },
                              getConditions = function(){
                                list(reference = .self$reference, condition = base::setdiff(.self$conditions, .self$reference) )
                              },
                              getResultTable = function(){
                                pvalues <- .self$getModPValuesCI()
                                #pvalues <- subset(pvalues, select = c("effectSize","p.ord","p.mod","q.ord","q.mod","log2FC"))

                                intensityWithNA <- merge(data.frame(nrNAs = .self$getNrNAs()),.self$proteinIntensity, by="row.names" )
                                rownames(intensityWithNA) <- intensityWithNA$Row.names
                                intensityWithNA$Row.names <- NULL

                                intensityWithNA <- merge(intensityWithNA, .self$getNormalized()$data, by="row.names", suffix = c(".raw", ".transformed"))
                                rownames(intensityWithNA) <- intensityWithNA$Row.names
                                intensityWithNA$Row.names <- NULL

                                intensityWithNA <- merge(pvalues,intensityWithNA, by="row.names")
                                rownames(intensityWithNA) <- intensityWithNA$Row.names
                                intensityWithNA$Row.names <-NULL

                                grpAvg <- .self$getNormalizedGrpAverages()

                                prots <- .self$proteinAnnotation
                                mm <- merge(grpAvg, intensityWithNA , by="row.names")
                                rownames(mm) <- mm$Row.names
                                mm$Row.names <-NULL
                                mm <- merge(prots, mm, by="row.names")
                                mm$Row.names <-NULL
                                return(mm)
                              },
                              getResultTableWithPseudo = function(){
                                "add pseudo p-values and pseudo fold changes"
                                ### compute pseudo p-values
                                results <- .self$getResultTable()
                                c2name <- base::setdiff(.self$conditions, .self$reference)


                                # fix references
                                r <-results[,.self$reference]
                                romit <- na.omit(r)
                                minr <- mean( sort(romit)[1:ceiling(length(romit)/10)] )

                                r[is.na(r)] <- minr
                                results[, paste("pseudo.", .self$reference, sep="")] <- r


                                # fix condition
                                c2 <- results[,c2name]
                                c2omit <- na.omit(c2)
                                minc2 <-  mean( sort(c2omit)[1:ceiling(length(c2omit)/10)])

                                c2[is.na(c2)] <- minc2
                                results[, paste("pseudo.", c2name, sep="")]<- c2


                                results$pseudo.log2FC <- c2 - r
                                results <- dplyr::mutate(results, pseudo.P.Value = ifelse(is.na(P.Value), 0, P.Value))
                                results <- dplyr::mutate(results, pseudo.adj.P.Val = ifelse(is.na(adj.P.Val), 0, adj.P.Val))

                                return(results)
                              },
                              getNormalizedGrpAverages = function(){
                                "computes grp averages per protein"
                                nrdata <-.self$getNormalizedConditionData(.self$getConditions()[1])
                                rowNas <- quantable::rowNAs(nrdata)
                                grpAverage1 <- apply(nrdata, 1, mean, na.rm = TRUE)
                                grpAverage1[rowNas == ncol(nrdata)] <- NA

                                nrdata <-.self$getNormalizedConditionData(.self$getConditions()[2])
                                rowNas <- quantable::rowNAs(nrdata)
                                grpAverage2 <- apply(nrdata, 1, mean, na.rm = TRUE)
                                grpAverage2[rowNas == ncol(nrdata)] <- NA

                                grpAverage<-cbind(grpAverage1,grpAverage2)
                                colnames(grpAverage) <- .self$getConditions()
                                return(grpAverage)
                              },
                              getResultTableWithPseudoAndClustering = function(){
                                "calls getResultTableWithPseudo() and adds clusterID column"

                                # Run getResultTableWithPseudo()
                                tmpdata <- .self$getResultTableWithPseudo()

                                # Retrieve normalized data for clustering
                                tmpdata1 <- .self$getNormalized()$data

                                # Clustering using quantable::simpleheatmap3()
                                clustering <- quantable::simpleheatmap3(tmpdata1,
                                                                        labCol = colnames(tmpdata1),
                                                                        labRow = rownames(tmpdata1),
                                                                        plot = FALSE,
                                                                        nrOfClustersRow = .self$numberOfProteinClusters)
                                row_clustering <- clustering$Row
                                colnames(row_clustering) = c("TopProteinName", "ClusterID")

                                # Merging ResultTable with ClusterID data.frame
                                results <- dplyr::full_join(tmpdata, row_clustering, by = "TopProteinName")
                                return(results)
                              },
                              getNumberOfClusters = function(){
                                .self$numberOfProteinClusters
                              }
                            )
)



