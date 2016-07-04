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
#' library(SRMService)
#' tmp <- "D:/googledrive/DataAnalysis/p1930"
#' allData <- read.csv( file=file.path(tmp,"data/longFormat.txt"),row.names = 1)
#' data <-allData
#' pool=2
#' data <-allData[allData$pool==pool,]
#' srms <- SRMService(data,qvalue=0.25)
#' srms$plotQValues()
#' srms$plotTransition()
#' srms$plotTransition(light=TRUE)
#' srms$getNrNAs()
#' srms$getNrNAs(light=TRUE)
#' srms$maxNAHeavy
#' srms$maxNALight
#' hist(srms$piw$nrNA)
#'
#' srms$setMaxNAHeavy(30)
#' srms$setMaxNALight(40)
#'
#' tmpH <- srms$getTransitionIntensities()
#' dim(tmpH)
#' stopifnot(max(apply(tmpH, 1, function(x){sum(is.na(x))}))<=30)
#' tmpL <- srms$getTransitionIntensities(light=TRUE)
#' dim(tmpL)
#' max(apply(tmpL, 1, function(x){sum(is.na(x))}))
#' x<-srms$getMatchingIntensities()
#' dim(x$light)
#' tmp<-srms$getLHLog2FoldChange()
#'
SRMService <- setRefClass("SRMService",
                          fields = list( data = "data.frame",
                                         dataq = "data.frame",
                                         qValueThreshold = "numeric",
                                         maxNAHeavy = "numeric",
                                         maxNALight = "numeric",
                                         piw = "data.frame",
                                         int = "data.frame"

                          ),methods = list(
                            setQ = function(qvalue = 0.05){
                              .self$dataq <- .self$data
                              .self$qValueThreshold <- qvalue
                              message("setting stuff")
                              .self$dataq$Area[.self$dataq$annotation_QValue > .self$qValueThreshold] <- NA
                              .makePivotData()
                            },
                            .makePivotData = function(){
                              message("pivoting data")
                              .self$piw <- SRMService::piwotPiw(.self$dataq)
                              .self$int <- SRMService::getIntensities(.self$piw)
                              rownames(.self$piw) <-rownames(.self$int)
                              .self$piw <- .self$piw[,1:6]

                              nas <- apply(.self$int , 1 , function(x){sum(is.na(x))})
                              .self$piw$nrNA <- nas

                            },
                            initialize = function(data,
                                                  qvalue = 0.05
                            ){
                              stopifnot(getRequiredColumns() %in% colnames(data))
                              .self$data <- data

                              .self$maxNAHeavy <- length(unique(.self$data$Replicate.Name))
                              .self$maxNALight <- length(unique(.self$data$Replicate.Name))
                              setQ(qvalue)
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
                              isolable<-ifelse(light, "light", "heavy")
                              nas <-subset(.self$piw,Isotope.Label==isolable)$nrNA
                              hist(nas, main=isolable)
                              abline(v=ifelse(light, .self$maxNALight, .self$maxNAHeavy), col=2)
                              invisible(nas)
                            },
                            setMaxNAHeavy = function(max=0){
                              "set maximum of na's heavy row"
                              .self$maxNAHeavy = max
                            },
                            setMaxNALight = function(max=0){
                              "set maximum of na's in light row"
                              .self$maxNALight = max
                            },
                            getTransitionIntensities=function(light=FALSE){
                              "get matrix with intensities, where nr of NAs in row < maxNA"
                              idx <- .self$piw$nrNA <= ifelse(light,.self$maxNALight, .self$maxNAHeavy)

                              int_ <- subset(.self$int, .self$piw$Isotope.Label==ifelse(light, "light", "heavy")
                                             & idx
                              )
                              return(int_)
                            },
                            getMatchingIntensities = function(){
                              "get matrix with intensities, where nr of NAs in row < maxNA"
                              intLight_ <- .self$getTransitionIntensities(light=TRUE)
                              intHeavy_ <- .self$getTransitionIntensities()

                              rownames(intLight_) <- gsub("_light$","",rownames(intLight_))
                              rownames(intHeavy_) <- gsub("_heavy$","",rownames(intHeavy_))

                              idx <-intersect(rownames(intLight_),rownames(intHeavy_))

                              return(list(light = intLight_[idx,], heavy = intHeavy_[idx,]))
                            },
                            plotCommonTransitions = function(light=FALSE){
                              "Shows transitions which occure in heavy and light"
                              int_ <- getMatchingIntensities()
                              int_ <- ifelse(light,int_$light, int_$heavy)
                              imageWithLabels(log2(t((int_))),col = quantable::getRedScale(),
                                              main="heavy Int",marLeft=c(5,15,3,3),marRight = c(5,0,3,3))
                              invisible(int_)
                            },
                            getLHLog2FoldChange = function(plot=TRUE){
                              int_<-getMatchingIntensities()
                              stopifnot(colnames( int_$light) == colnames(int_$heavy))
                              stopifnot(rownames( int_$light) == rownames(int_$heavy))

                              logfc <-(log2(( int_$light )) - log2((int_$heavy)) )

                              if(plot){
                                imageWithLabels(t(logfc), main="log2(L/H)", col= getBlueWhiteRed(), marLeft=c(5,15,3,3),marRight = c(5,0,3,3))
                                invisible(logfc)
                              }
                              else{
                                return(logfc)
                              }
                            },

                            plotTransition = function(light = FALSE){
                              int_ <- .self$getTransitionIntensities(light=light)
                              quantable::imageWithLabels( t(as.matrix(log2( int_ ) )) , col = quantable::getRedScale(),
                                                          main=ifelse(light, "ligth", "heavy"),marLeft=c(5,15,3,3),marRight = c(5,0,3,3))
                            }

                          )
)



