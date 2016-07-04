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
#' pool=1
#' data <-allData[allData$pool==pool,]
#' colnames(data)
#' length(unique(data$Replicate.Name))
#' srms <- SRMService(data,qvalue=0.25)
#' srms$plotQValues()
#' srms$plotTransition()
#' srms$plotTransition(light=TRUE)
#' srms$getNrNAs()
#' srms$getNrNAs(light=TRUE)
#' tmpH <- srms$getTransitionIntensities()
#' head(tmpH)
#' tmpL <- srms$getTransitionIntensities(light=TRUE)
#'
#' head(tmpL)
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

                            getTransitionIntensities=function(light=FALSE){
                              int_ <- subset(.self$int, .self$piw$Isotope.Label==ifelse(light, "light", "heavy"))
                              return(int_)
                            },

                            plotTransition = function(light = FALSE){
                              int_ <- .self$getTransitionIntensities(light=light)
                              quantable::imageWithLabels( t(as.matrix(log2( int_ ) )) , col = quantable::getRedScale(),
                                                          main=ifelse(light, "ligth", "heavy"),marLeft=c(5,15,3,3),marRight = c(5,0,3,3))
                            }

                          )
)



