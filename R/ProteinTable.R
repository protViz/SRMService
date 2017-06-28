#' Protein Table
#' @docType class
#' @export
#' @field data data.frame with colnames sample ID rownames proteinID
#' @field conditionmapping data.frame with 2 columns providing mapping of sampleID to condition
#' @field experimentID name of the experiment
#' @importFrom R6 R6Class
#'
ProteinTableR6 <- R6::R6Class("ProteinTableR6",
                        public = list(
                          data = NULL,
                          conditionmap = NULL,
                          experimentID = NULL,
                          housekeeper = NULL,
                          normalized=NULL,

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

                          condition_colors=function(){
                            fact <- as.factor(self$conditionmap[colnames(self$data), "Condition"])
                            tmpcol <- as.numeric(fact)
                          },
                          setHouseKeepers = function(housekeeper){
                            self$housekeeper <- housekeeper[housekeeper %in% rownames(self$data)]
                            diff <- setdiff(housekeeper,self$housekeeper)
                            if( length(diff) > 0){
                              warning(diff, " are not among the proteins")
                            }
                            invisible(self$housekeeper)
                          },
                          normalizeUsingHouseKeeping = function( protein , FUN = median , plot=TRUE ){
                            "normalize using housekeeping proteins"
                            if(missing(protein)) {
                              .housekeepers <- self$data[self$housekeeper,]
                              normalize <- apply(.housekeepers,2, FUN, na.rm=T)
                            } else {
                              .housekeepers <- self$data[protein,]
                              normalize <- unlist(.housekeepers)
                            }
                            if(sum(is.na(normalize))>0){
                              warning("can not normalize data, protein not quantified in some samples. NAs")
                              return(FALSE)
                            }
                            print(rownames(.housekeepers))

                            if(plot){
                              nrlines <- nrow(.housekeepers)
                              matplot(t(.housekeepers),lwd=2,
                                      col=1:nrlines,
                                      lty=1:nrlines,
                                      type="l",las=2,xaxt="n")
                              axis(1,at=1:ncol(.housekeepers), labels=colnames(.housekeepers),las=2)
                              legend("bottomleft",legend=(colnames(t(.housekeepers))),col=1:nrlines, lty=1:nrlines , lwd=2)
                              lines(1:length(normalize), normalize, col=2)
                            }

                            normprotquant <- sweep(self$data, 2, normalize, "-" )
                            pp <-ProteinTable(normprotquant,
                                              self$conditionmap,
                                              experimentID = self$experimentID
                            )
                            pp$setHouseKeepers(rownames(.housekeepers))
                            pp$normalized = rownames(.housekeepers)
                            return(pp)
                          },
                          getProtein=function(protein){
                            return(self$data[protein,])
                          },
                          getProteinMatrix=function(){
                            res <- self$data
                            colnames(res) <-  self$conditionmap[colnames(self$data),"Colnames"]
                            invisible(res)
                          }
                        )
)
