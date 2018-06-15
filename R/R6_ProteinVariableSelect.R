library(R6)

#' Wraps 4 column matrix, with Condition, Run, Protein, Intensity
#' @importFrom reshape2 dcast
#' @export
ProteinVariableSelect <- R6Class("ProteinVariableSelect",
                                 public = list(
                                   data = NULL,
                                   results = list(),
                                   initialize = function(data){
                                     stopifnot(c("Condition", "Run", "Protein", "Intensity") %in% colnames(data))
                                     self$data <- data
                                   }
                                   ,getIntensities = function(){
                                     datamI <- self$getWideFormat()
                                     datamI <- datamI[,3:ncol(datamI)]
                                     rownames(datamI) <-  self$getWideFormat()$Run
                                     return(datamI)
                                   }
                                   ,getWideFormat = function(){
                                     res <- dcast(self$data, Run + Condition ~ Protein, value.var = "Intensity")
                                     res
                                   }
                                   ,getWideNoMissing = function(){
                                     intensities <- self$getIntensities()
                                     medians <- apply( intensities , 2 , median, na.rm = TRUE )
                                     for(i in 1:ncol(intensities)){
                                       intensities[is.na(intensities[,i]),i]<- medians[i]
                                     }
                                     return(intensities)

                                   },
                                   getIntensitiesNoMissing=function(){
                                     self$getWideNoMissing()
                                   }
                                 ))







