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
#' head(allData)
#' data <-allData
#' pool=2
#' data <-allData[allData$pool==pool,]
#' head(data)
#' srms <- SRMService(data,qvalue=0.05)
#' head(srms$piw)
#'
#' srms$plotQValues()
#' srms$plotTransition()
#' srms$plotTransition(light=TRUE)
#'
#' srms$getNrNAs()
#' srms$getNrNAs(light=TRUE)
#' srms$maxNAHeavy
#' srms$maxNALight
#'
#' srms$plotCommonTransitions()
#' srms$getLHLog2FoldChange()
#' srms$setMaxNAHeavy(30)
#' srms$setMaxNALight(40)
#'
#' colnames(srms$piw)
#' tmpH <- srms$getTransitionIntensities()
#' dim(tmpH)
#' stopifnot(max(apply(tmpH, 1, function(x){sum(is.na(x))}))<=30)
#' tmpL <- srms$getTransitionIntensities(light=TRUE)
#' dim(tmpL)
#' stopifnot(max(apply(tmpL, 1, function(x){sum(is.na(x))}))<=40)
#' x<-srms$getMatchingIntensities()
#' names(x)
#' dim(x$light)
#' dim(x$heavy)
#'
SRMService <- setRefClass("SRMService",
                          fields = list( data = "data.frame",
                                         dataq = "data.frame",
                                         qValueThreshold = "numeric",
                                         maxNAHeavy = "numeric",
                                         maxNALight = "numeric",
                                         piw = "data.frame",
                                         int = "data.frame",
                                         lightLabel = "character",
                                         heavyLable = "character"


                          ),methods = list(
                            initialize = function(data,
                                                  qvalue = 0.05
                            ){
                              .self$lightLabel = "light"
                              .self$heavyLable = "heavy"

                              stopifnot(getRequiredColumns() %in% colnames(data))
                              .self$data <- data[,getRequiredColumns()]
                              .self$maxNAHeavy <- length(unique(.self$data$Replicate.Name))
                              .self$maxNALight <- length(unique(.self$data$Replicate.Name))
                              setQ(qvalue)
                            },
                            setQ = function(qvalue = 0.05){
                              .self$dataq <- .self$data
                              .self$qValueThreshold <- qvalue
                              message( "Setting intensities to NA for qvalues larger than: ", .self$qValueThreshold )
                              .self$dataq$Area[.self$dataq$annotation_QValue > .self$qValueThreshold] <- NA
                              .makePivotData()
                            },
                            .mergeHL=function(piwdata){
                              " make sure that to every light you have also an heavy transition "
                              d2 <- reshape2::melt(piwdata, id.vars= colnames(piwdata)[1:6], variable.name = 'Replicate.Name',value.name='Area')
                              dl <- d2[d2$Isotope.Label == .self$lightLabel,]
                              dh <- d2[d2$Isotope.Label == .self$heavyLable,]

                              dl <- subset(dl, select = colnames(dl)!="Isotope.Label")
                              dh <- subset(dh, select = colnames(dh)!="Isotope.Label")

                              reportMissing(dl,dh)

                              colnames(dl)[colnames(dl) == "Area"] <- "light"
                              colnames(dh)[colnames(dh) == "Area"] <- "heavy"


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
                              isolable <- ifelse(light, "light", "heavy")
                              nas <- subset(.self$piw,Isotope.Label==isolable)$nrNA
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
                              int_ <- if(light){ int_$light}else{int_$heavy}
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



