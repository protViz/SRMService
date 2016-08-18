library(quantable)
.reportMissing <- function(dl, dh){
  df_args <- c( subset(dl, select = colnames(dl)!="Area"), sep=".")
  dlid <-do.call(paste, df_args)
  df_args <- c( subset(dh, select = colnames(dh)!="Area"), sep=".")
  dhid <-do.call(paste, df_args)

  missingInHeavy <- setdiff(dlid,dhid)
  missingInLight <- setdiff(dhid,dlid)
  if(length(missingInHeavy) > 0){
    warning("Transitions present in Light but missing in Heavy")
    warning(missingInHeavy)
  }
  if(length(missingInLight) > 0){
    warning("Transitions present in Light but missing in Light")
    warning(missingInHeavy)
  }
}

.findDecorrelated <- function(res, threshold = 0.8){
  nrtrans <- ncol(res)
  ids <- rowSums(res < threshold)
  names(which((nrtrans-1)== ids))
}

#' Peptide table
#' @export PeptideTable
#' @exportClass PeptideTable
PeptideTable <- setRefClass("PeptideTable",
                           fields = list( data = "data.frame",
                                          ids = "data.frame",
                                          experimentID= "character")
                           ,methods = list(
                             initialize = function(data, ids, experimentID){
                               .self$experimentID = experimentID
                               .self$data = data
                               .self$ids = ids
                               rownames(.self$ids) <- rownames(.self$data)
                             },
                             getProteinsAsList = function(){
                               xx <- .self$ids$Protein.Name
                               peptab <- by(.self$data ,INDICES=xx,function(x){x})
                               return(peptab)
                             },
                             getProteinIntensities = function(plot=TRUE, FUN = median, scale=TRUE){
                               proteins <- (aggregate(.self$data, list(Protein.Name=.self$ids$Protein.Name),FUN, na.rm=TRUE))
                               rownames(proteins) <- proteins$Protein.Name
                               proteins <- proteins[,2:ncol(proteins)]
                               if(plot){
                                 toplot <- if(scale){scale(t(proteins))}else{t(proteins)}
                                 imageWithLabels(toplot ,
                                                 main="log2(L/H)",
                                                 col= getBlueWhiteRed(),
                                                 marLeft=c(5,15,3,3),
                                                 marRight = c(5,0,3,3))
                               }
                               invisible(proteins)
                             },
                             plot = function(){
                                 imageWithLabels(t(.self$data) ,
                                                 main="log2(L/H)",
                                                 col= getBlueWhiteRed(),
                                                 marLeft=c(5,15,3,3),
                                                 marRight = c(5,0,3,3))

                             }
                           )
)

#' Transition table
#'
#' @export TransitionTable
#' @exportClass TransitionTable
TransitionTable <- setRefClass("TransitionTable",
                               fields = list( data = "data.frame",
                                              ids = "data.frame",
                                              experimentID = "character")
                               ,methods = list(
                                 initialize = function(data, experimentID=""){
                                   .self$experimentID = experimentID
                                   .self$data = data
                                   .self$ids <- .self$rownamesAsTable(rownames(data))
                                   rownames(.self$ids) <- rownames(.self$data)
                                 },
                                 rownamesAsTable = function(x){
                                   "help function"
                                   tab <- data.frame(quantable::split2table(x,split="\\_"))
                                   colnames(tab) <- getIDLabels()[1:ncol(tab)]
                                   return(tab)
                                 },
                                 filterData = function(minNrTransition = 2, minNrPeptides = 0){
                                   trans <- aggregate(rep(1, nrow(.self$ids)),
                                                      by = list(Peptide.Sequence = .self$ids$Peptide.Sequence,
                                                                Precursor.Charge = .self$ids$Precursor.Charge ),
                                                      length)
                                   trans <- trans[trans$x >= minNrTransition,]
                                   .ids <- merge(trans[,1:2], .self$ids)
                                   .ids <- .ids[,getIDLabels()[1:ncol(.ids)]]
                                   df_args <- c(.ids, sep="_")
                                   .ids <- do.call(paste, df_args)
                                   return(TransitionTable(.self$data[.ids, ], .self$experimentID))
                                 },
                                 dim=function(){
                                   dd <- base::dim(.self$data)
                                   return(dd)
                                 },
                                 getCorrelatedPeptides = function(is = TRUE, minCorrelation = 0.8){
                                   ' returns peptides with well correlated transitions'
                                   prottab <- .self$getPeptidesAsList()
                                   res<-vector(mode = "list", length(prottab))
                                   for(i in 1:length(prottab)){

                                     res[[i]] <- SRMService::transitionCorrelations(dataX)
                                   }
                                   if(is){
                                     idx <-which(sapply(res, min) >= 0.8)
                                   }else{
                                     idx <-which(sapply(res, min) < 0.8)
                                   }
                                   return(prottab[idx])
                                 }, getPeptidesAsList = function(){
                                   ' returns list of matrices each matrix representing single peptide'

                                   xx <- .self$getPeptideIDs()
                                   paste_args <- c(xx, sep="_")
                                   prottab <- by(.self$data ,INDICES=do.call(paste,paste_args),function(x){x})
                                   return(prottab)
                                 },
                                 removeDecorrelated = function(minCorrelation = 0.8){
                                   "removes decorrelated peptides"
                                   prottab <- .self$getPeptidesAsList()
                                   res <- lapply(prottab, SRMService::transitionCorrelations)
                                   toremove <- unlist(lapply(res, .findDecorrelated, threshold=minCorrelation))
                                   xx<-TransitionTable(filteredfc$data[setdiff(rownames(filteredfc$data),toremove),],.self$experimentID)
                                   return(xx)
                                 },
                                 plot = function(){
                                   imageWithLabels(.self$data , main="log2(L/H)",
                                                   col= getBlueWhiteRed(),
                                                   marLeft=c(5,15,3,3),
                                                   marRight = c(5,0,3,3))
                                 },
                                 getPeptideIDs = function(){
                                   .self$ids[,c("Protein.Name","Peptide.Sequence","Precursor.Charge")]
                                 },
                                 getPeptideIntensities = function(plot=TRUE, FUN = median){
                                   peptides <- (aggregate(.self$data, .self$getPeptideIDs(),FUN, na.rm=TRUE))
                                   rownames(peptides) <- do.call(paste, c(peptides[,1:ncol(.self$getPeptideIDs())], sep="_"))


                                   invisible(PeptideTable(data=peptides[,(ncol(.self$getPeptideIDs())+1):ncol(peptides)],
                                                          ids = peptides[,1:ncol(.self$getPeptideIDs())]), .self$experimentID)
                                 },
                                 getProteinsAsList = function(){
                                   xx <- .self$ids$Protein.Name
                                   peptab <- by(.self$data ,INDICES=xx,function(x){x})
                                   return(peptab)
                                 }


                               ))
#' R access to Bibliospec File
#'
#' @description
#' This class implements an R referenz class for BiblioSpec
#' generated sqlite files and can return the data contained as data.frames
#' or as list of tandem mass spectra peptide assignments objects (psm).
#'
#'
#' @field dbfile database file location
#' @import parallel
#' @import DBI
#' @import quantable
#' @import methods
#' @export SRMService
#' @exportClass SRMService
#' @details The function performs a SQL query on the SQLite
#'
#' @examples
#'
#' library(SRMService)
#' tmp <- "D:/googledrive/DataAnalysis/p1930"
#' allData <- read.csv( file=file.path(tmp,"data/longFormat.txt"),row.names = 1)
#' head(allData)
#' data <-allData
#' pool=2
#' data <-allData[allData$pool==pool,]
#' head(data)
#' srms <- SRMService(data,qvalue=0.05)
#' SRMService$methods()
#' srms$maxNAHeavy()
#' SRMService$fields()
#' head(srms$piw)
#'
#' srms$plotQValues()
#' srms$plotTransitions()
#' srms$plotTransitions(light=TRUE)
#'
#' srms$getNrNAs()
#' srms$getNrNAs(light=TRUE)
#' srms$maxNAHeavy()
#' srms$maxNALight()
#' srms$maxNAHeavy(80)
#' srms$maxNALight(80)
#' srms$plotCommonTransitions()
#' srms$plotCommonTransitions(light=TRUE)
#' tmp <-srms$getLHLog2FoldChange(maxNA = 20)
#' srms$maxNAFC()
#' resH <- srms$plotCommonTransitions()
#' dim(resH)
#' resL <-srms$plotCommonTransitions(light=TRUE)
#' dim(resL)
#' resAll <- srms$getMatchingIntensities()
#' dim(resAll$light)
#'
#' tmpH <- srms$getTransitionIntensities()$data
#' tmpL <- srms$getTransitionIntensities(light=TRUE)$data
#'
#' dim(tmpL)
#'
SRMService <- setRefClass("SRMService",
                          fields = list( data = "data.frame",
                                         dataq = "data.frame",
                                         qValueThreshold = "numeric",
                                         MaxNAHeavy = "numeric",
                                         MaxNALight = "numeric",
                                         MaxNAFC = "numeric",
                                         piw = "data.frame",
                                         int = "data.frame",
                                         lightLabel = "character",
                                         heavyLabel = "character",
                                         experimentID = "character"
                          ),methods = list(
                            initialize = function(data,experimentID="",
                                                  qvalue = 0.05
                            ){
                              .self$experimentID = experimentID
                              .self$lightLabel = "light"
                              .self$heavyLabel = "heavy"

                              stopifnot(getRequiredColumns() %in% colnames(data))
                              .self$data <- data[,getRequiredColumns()]
                              .self$MaxNAHeavy <- length(unique(.self$data$Replicate.Name))
                              .self$MaxNALight <- .self$MaxNAHeavy
                              .self$MaxNAFC <- .self$MaxNAHeavy
                              setQ(qvalue)
                            },
                            setQ = function(qvalue = 0.05){
                              .self$dataq <- .self$data
                              .self$qValueThreshold <- qvalue
                              message( "Setting intensities to NA for qvalues larger than: ", .self$qValueThreshold )
                              .self$dataq$Area[.self$dataq$annotation_QValue > .self$qValueThreshold] <- NA
                              .makePivotData()
                            },
                            qValueHist=function(){
                              hist(.self$data$annotation_QValue, main="q Values")
                              abline(v = .self$qValueThreshold, col=2)
                            },
                            .mergeHL=function(piwdata){
                              library(reshape2)
                              " make sure that to every light you have also an heavy transition "
                              d2 <- reshape2::melt(piwdata, id.vars= colnames(piwdata)[1:6], variable.name = 'Replicate.Name',value.name='Area')
                              dl <- d2[d2$Isotope.Label == .self$lightLabel,]
                              dh <- d2[d2$Isotope.Label == .self$heavyLabel,]

                              dl <- subset(dl, select = colnames(dl)!="Isotope.Label")
                              dh <- subset(dh, select = colnames(dh)!="Isotope.Label")

                              .reportMissing(dl,dh)

                              colnames(dl)[colnames(dl) == "Area"] <- .self$lightLabel
                              colnames(dh)[colnames(dh) == "Area"] <- .self$heavyLabel


                              fixedData <- merge(dl,dh)
                              tmp <-melt(fixedData, id.vars = colnames(fixedData)[1:6],variable.name = "Isotope.Label",value.name = "Area" )
                              return(tmp)
                            }
                            ,.makePivotData = function(){
                              message("pivoting data")
                              .self$piw <- SRMService::piwotPiw(.self$dataq)
                              .self$dataq <- .mergeHL(.self$piw)
                              print("mergeDone!")
                              .self$piw <- SRMService::piwotPiw(.self$dataq)
                              .self$int <- SRMService::getIntensities(.self$piw)
                              rownames(.self$piw) <-rownames(.self$int)
                              nas <- apply(.self$int , 1 , function(x){sum(is.na(x))})
                              .self$piw$nrNA <- nas
                            },

                            summary = function(){
                              "summarize experiment"
                            },
                            plotQValues  = function(){
                              "show q value distribution for peak groups"
                              hist(.self$data$annotation_QValue)
                              abline(v=.self$qValueThreshold,col=2)
                            },
                            getNrNAs = function(light = FALSE){
                              "show nr of NA's for heavy (defalt) or light transitions"
                              isolabel <- if(light){ .self$lightLabel} else {.self$heavyLabel}
                              nas <- subset(.self$piw,Isotope.Label==isolabel)$nrNA
                              plot(table(nas), main=isolabel)
                              abline(v=ifelse(light, .self$MaxNALight, .self$MaxNAHeavy), col=2)
                              invisible(nas)
                            },
                            maxNAHeavy = function(max){
                              "set maximum of na's heavy row"
                              if(!missing(max)){
                                .self$MaxNAHeavy = max
                              }else{
                                return(.self$MaxNAHeavy)
                              }
                            },
                            maxNALight = function(max){
                              "set maximum of na's in light row"
                              if(!missing(max)){
                                .self$MaxNALight = max
                              }else{
                                return(.self$MaxNALight)
                              }
                            },
                            maxNAFC = function(max){
                              if(!missing(max)){
                                .self$MaxNAFC = max
                              }else{
                                return(.self$MaxNAFC)
                              }
                            },
                            getTransitionIntensities=function(maxNA,light=FALSE){
                              "get matrix with intensities, where nr of NAs in row < maxNA"
                              if(!missing(maxNA)){
                                idx <- .self$piw$nrNA <= maxNA
                              }else{
                                idx <- .self$piw$nrNA <= ifelse(light,.self$MaxNALight, .self$MaxNAHeavy)
                              }
                              int_ <- subset(.self$int,
                                             .self$piw$Isotope.Label==ifelse(light,.self$lightLabel, .self$heavyLabel)
                                             & idx
                              )
                              return(TransitionTable(int_,.self$experimentID))
                            }
                            ,
                            stripLabel=function(int_,light=FALSE){
                              label<-if(light){.self$lightLabel}else{.self$heavyLabel}
                              rownames(int_) <- gsub(paste("_",label,"$",sep=""),"",rownames(int_))
                              invisible(int_)
                            }
                            ,
                            getMatchingIntensities = function(){
                              "get matrix with intensities, where nr of NAs in row < maxNA"
                              intLight_ <- .self$getTransitionIntensities(light=TRUE)$data
                              intHeavy_ <- .self$getTransitionIntensities()$data

                              intLight_ <- stripLabel(intLight_, light=TRUE)
                              intHeavy_ <- stripLabel(intHeavy_)

                              idx <- intersect(rownames(intLight_),rownames(intHeavy_))

                              return(list(light = intLight_[idx,], heavy = intHeavy_[idx,]))
                            },
                            plotCommonTransitions = function(light=FALSE){
                              "Shows transitions which occure in heavy and light"
                              int_ <- getMatchingIntensities()
                              int_ <- if(light){ int_$light}else{int_$heavy}
                              main <- paste(if(light){.self$lightLabel} else {.self$heavyLabel}, "Intensity")
                              imageWithLabels(log2(t((int_))),col = quantable::getRedScale(),
                                              main=main,marLeft=c(5,15,3,3),marRight = c(5,0,3,3))
                              invisible(int_)
                            },
                            getLHLog2FoldChange = function( maxNA, minNrOfTransitons=2, plot=TRUE){
                              int_<-getMatchingIntensities()
                              stopifnot(colnames( int_$light) == colnames(int_$heavy))
                              stopifnot(rownames( int_$light) == rownames(int_$heavy))

                              logfc <-(log2(( int_$light )) - log2((int_$heavy)) )
                              if(missing(maxNA)){
                                maxNA <- .self$MaxNAFC
                              }
                              logfc <- subset(logfc , maxNA >= apply(logfc ,1, function(x){sum(is.na(x))}))

                              if(plot){
                                imageWithLabels(t(logfc), main=
                                                  paste("log2(",.self$lightLabel, "/", .self$heavyLabel, ")"),
                                                col= getBlueWhiteRed(),
                                                marLeft=c(5,15,3,3),
                                                marRight = c(5,0,3,3))
                              }
                              invisible(TransitionTable(logfc,.self$experimentID))

                            },

                            plotTransitions = function(light = FALSE){
                              int_ <- .self$getTransitionIntensities(light=light)$data
                              quantable::imageWithLabels( t(as.matrix(log2( int_ ) )) ,
                                                          col = quantable::getRedScale(),
                                                          main=ifelse(light, .self$lightLabel,.self$heavyLabel),
                                                          marLeft=c(5,15,3,3),
                                                          marRight = c(5,0,3,3))
                            }

                          )
)



