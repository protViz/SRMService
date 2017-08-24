#' Condition Map
#' @docType class
#' @export
#' @include RequiredColumns.R
#' @examples
#'
#' x <- data.frame(Raw.file = c("A","B","C"),"MeasOrder" = c(1,2,3) , Gender = c("F","F","M"))
#' cm <- Annotation$new(x, fixed = "Gender")
#' cm$exists("Gender")
#' cm$exists("Stop")
#' cm$fixed
#' cm$random
#' cm$as_numeric("Gender")
#' cm$levels("Gender")
#' cm$get("Gender")
#' cm$get_color("Gender")
#' cm$get_required()
Annotation <- R6::R6Class("Annotation",
                          public = list(
                            annotation = data.frame(),
                            fixed = NULL,
                            random = NULL,
                            required = c("Raw.file","MeasOrder"),
                            initialize = function(annotation, fixed =NULL , random =NULL){
                              if(sum(self$required %in% colnames(annotation))!=length(self$required)){
                                stop("condition mappings does not contain columns : ", self$required)
                              }
                              self$fixed = c(fixed)
                              self$random = c(random)
                              if( !is.null(self$fixed) ){
                                if( length(setdiff(self$fixed, colnames(annotation))) != 0 ){
                                  stop("annotation does not contain fixed effects : ", setdiff(self$fixed, colnames(annotation)))
                                }
                              }
                              if( !is.null(self$random)){
                                if( length(setdiff(self$random, colnames(annotation))) != 0 ){
                                  stop("annotation does not contain fixed effects : ", setdiff(self$random, colnames(annotation)))
                                }
                              }
                              self$annotation = annotation
                            },
                            exists = function(colname){
                              if(length(setdiff(colname , colnames(self$annotation)))){
                                stop("annotation does not contain  : ", setdiff(colname, colnames(self$annotation)))
                              }
                              return(TRUE)
                            },
                            levels = function(colname){
                              self$exists(colname)
                              if(length(colname) > 1){
                                stop("can get levels only for a singl column")
                              }
                              levels(as.factor(unlist(self$get(colname ))))
                            },
                            get = function(colname){
                              self$exists(colname)
                              subset( self$annotation, select = colname )
                            },
                            get_factors = function(){
                              return(colnames(self$annotation))
                            },
                            get_required = function(){
                              return(subset(self$annotation, select = self$required))
                            },
                            as_numeric = function(colname){
                              self$exists(colname)
                              as.numeric(as.factor( unlist(self$get(colname)) ))
                            },
                            get_samples = function(){
                              subset(self$annotation, select = self$required[1])
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
#' @field conditionmapping data.frame with 2 columns providing mapping of sampleID to condition
#' @field experimentID name of the experiment
#' @importFrom R6 R6Class
#' @examples
#' pt <- ProteinTableR6$new()
ProteinTableR6 <- R6::R6Class("ProteinTableR6",public = list(
  proteinIntensities = NULL,
  required = c("ProteinID", "NrPeptides"),
  annotation = NULL,
  experimentID = NULL,

  initialize = function(data, annotation, experimentID="" ){
    self$annotation = annotation

    check <- setdiff(colnames(data) , annotation$get_samples())
    if(length(check)!=0){
      warning(check)
      stop("Colnames in data do not match sampleID")
    }
    self$experimentID = experimentID
    self$proteinIntensities = data

  },
  getProtein=function(protein){
    "return data for protein"
    return(subset(self$data, proteinID = protein))
  },
  getLongFormat = function(){

  },
  get_annotation = function(){
    self$annotation
  }
)
)
