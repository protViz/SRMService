#' Condition Map
#' @docType class
#' @export
#' @import RequiredColumns.R
#'
ConditionMap <- R6::R6Class("ConditionMap",public = list(
  conditionmap = data.frame(),
  required = SRMService::annotationColumns,
  condition_colors=function(){
    fact <- as.factor(self$conditionmap[colnames(self$data), "Condition"])
    tmpcol <- as.numeric(fact)
  }

)
)


#' Protein Table
#' @docType class
#' @export
#' @field data data.frame with colnames sample ID rownames proteinID
#' @field conditionmapping data.frame with 2 columns providing mapping of sampleID to condition
#' @field experimentID name of the experiment
#' @importFrom R6 R6Class
#' @examples
#' pt <- ProteinTableR6$new()
ProteinTableR6 <- R6::R6Class("ProteinTableR6",
                        public = list(
                          proteinIntensities = matrix(),
                          conditionmap = NULL,
                          experimentID = character(),

                          initialize = function(data, conditionmapping, experimentID="" ){
                            self$normalized = ""
                            reccolumns <- c("Condition","BioReplicte", "Replicate.Name")
                            if(sum(reccolumns %in% colnames(conditionmapping))!=2){
                              stop("condition mappings does not contain columns : ", reccolumns)
                            }
                            self$conditionmap = .fixConditionMapping(conditionmapping)

                            check <- setdiff(colnames(data) , conditionmapping$Replicate.Name)
                            if(length(check)!=0){
                              warning(check)
                              stop("Colnames data do not match conditionmappings")
                            }
                            self$experimentID = experimentID
                            self$data = data

                          },
                          getProtein=function(protein){
                            "return data for protein"
                            return(self$data[protein,])
                          },
                          getProteinMatrix=function(){
                            res <- self$data
                            colnames(res) <-  self$conditionmap[colnames(self$data),"Colnames"]
                            invisible(res)
                          }
                          getLongFormat = function(){

                          }
                          getConditionMap = function(){
                            self$conditionmap
                          }
                        )
)
