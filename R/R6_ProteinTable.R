#' Condition Map
#' @docType class
#' @export
#' @include RequiredColumns.R
#' @import scales
#' @importFrom tibble tibble
#' @import dplyr
#' @examples
#'
#' library(SRMService)
#' library(dplyr)
#' x <- data.frame(Raw.file = c("A","B","C"),"MeasurementOrder" = c(1,2,3) , Gender = c("F","F","M"))
#' cm <- Annotation$new(x, experimentID = "E1 with random Effect", fixed = "Gender", random="MeasurementOrder")
#' cm$exists("Gender")
#' cm$exists("Stop")
#' cm$fixed
#' cm$random
#' cm$as_numeric("Gender")
#' cm$levels("Gender")
#' cm$get("Gender")
#' cm$get_color("Gender")
#' cm$sampleID
#' cm$get_sample_id()
#' cm$annotation
Annotation <- R6::R6Class("Annotation",
                          public = list(
                            annotation = tibble::tibble(),
                            fixed = NULL,
                            random = NULL,
                            sampleID = NULL,
                            experimentID = NULL,
                            initialize = function(annotation,
                                                  experimentID,
                                                  fixed =NULL,
                                                  random =NULL,
                                                  sampleID = "Raw.file")
                            {
                              self$sampleID = sampleID
                              self$experimentID = experimentID
                              self$fixed = c(fixed)
                              self$random = c(random)
                              if(length(setdiff(c(self$sampleID,self$fixed,self$random), colnames(data))) == 0  ){
                                stop("condition mappings does not contain columns : ", c(self$sampleID,self$fixed,self$random), sep = "\n")
                              }
                              if( !is.null(self$fixed) ){
                                if( length(setdiff(self$fixed, colnames(annotation))) != 0 ){
                                  stop("annotation does not contain fixed effects : ", setdiff(self$fixed, colnames(annotation)))
                                }
                              }
                              if( !is.null(self$random)){
                                if( length(setdiff(self$random, colnames(annotation))) != 0 ){
                                  stop("annotation does not contain random effects : ", setdiff(self$random, colnames(annotation)))
                                }
                              }
                              self$annotation <- annotation %>%  mutate_if(is.factor, as.character)
                            },
                            exists = function(colname)
                            {
                              if(length(setdiff(colname , colnames(self$annotation)))){
                                warning("annotation does not contain  : ", setdiff(colname, colnames(self$annotation)))
                                return(FALSE)
                              }
                              return(TRUE)
                            },
                            levels = function(colname){
                              if(self$exists(colname)){
                                if(length(colname) > 1){
                                  stop("can get levels only for a single column")
                                }
                                levels(as.factor(unlist(self$get(colname ))))
                              }
                            },
                            get = function(colname){
                              if(self$exists(colname)){
                                subset( self$annotation, select = colname )
                              }
                            },
                            get_factors = function(){
                              return(colnames(self$annotation))
                            },
                            get_sample_id = function(){
                              return(as.character(unlist(subset(self$annotation, select = self$sampleID))))
                            },
                            as_numeric = function(colname){
                              if(self$exists(colname)){
                                as.numeric(as.factor( unlist(self$get(colname)) ))
                              }
                            },
                            get_color = function( colname, pal_function = dichromat_pal, pal_name= "BrowntoBlue.10" ){
                              "convenience wrapper to scales"
                              self$exists(colname)
                              n <- length(self$levels(colname))
                              colors <- pal_function(pal_name)(n)[self$as_numeric(colname)]
                              return(colors)
                            }
                          ))

#' Protein Table
#' @docType class
#' @export
#' @field data data.frame with colnames sample ID rownames proteinID
#' @field experimentID name of the experiment
#' @importFrom R6 R6Class
#' @import dplyr tidyr rlang
#' @examples
#' library(SRMService)
#' x <- data.frame(Raw.file = c("A","B","C"),"MeasurementOrder" = c(1,2,3) , Gender = c("F","F","M"))
#' cm <- Annotation$new(x, experimentID = "E1 fixed Gender", fixed = "Gender")
#' dataProt <- data.frame("proteinID" = c("P1","P2","P3"), NrPeptides = c(1,2,3), "Fasta.Headers" = rep("",3), "A" = rexp(3), "B" = rexp(3), "C" = rexp(3))
#' length(intersect("Annotation",class(dataProt))) > 0
#' cm$sampleID
#' pt <- ProteinTableR6$new(cm)
#' pt$set_data_wide(dataProt)
#' pt$get_data_long()
#' pt$get_data_wide()
#' pt$get_all_long()
ProteinTableR6 <- R6::R6Class( "ProteinTableR6", public = list(
  data = data.frame(),
  proteinID = NULL,
  annotation = NULL,
  required = NULL,

  initialize = function(annotation,
                        proteinID = "proteinID",
                        required = c( "NrPeptides", "Fasta.Headers")){
    if(!length(intersect("Annotation",class(annotation))) > 0)
    {
      stop("annotation parameter must be of class Annotation")
    }
    self$annotation = annotation
    self$required = required
    self$proteinID = proteinID
  },
  set_data_wide = function(data){
    all_required <- c(self$proteinID, self$required, self$annotation$get_sample_id())
    check <- setdiff(all_required , colnames(data))
    if(length(check) != 0){
      stop("not having required columns : " , paste(check , sep=" ")  )
    }
    # take minimal set of columns
    data <- dplyr::select(data, all_required )
    # convert to wide format
    self$data <- tidyr::gather(data, key=rlang::UQ(rlang::sym(self$annotation$sampleID)), value="Intensity", self$annotation$get_sample_id() )
  },
  get_protein_id = function(){
    dplyr::select(self$data, self$proteinID)
  },
  get_data_long = function(){
    self$data
  },
  get_all_long = function(){
    dplyr::inner_join(self$annotation$annotation, self$data, by = self$annotation$sampleID)
  },
  get_data_wide = function(){
    tidyr::spread(self$data, self$annotation$sampleID , "Intensity")
  },
  get_annotation = function(){
    self$annotation
  }
)
)
