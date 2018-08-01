# function for modelling go here.
#' rocs helper function
rocs <- function(data ,response, predictor){
  responseX <- data %>% pull(!!sym(response))
  predictorX <- data %>% pull(!!sym(predictor))
  levels = levels(as.factor(responseX))
  if(length(levels) < 2){
    return(NULL)
  }
  res <- list()
  comparisons <- combn(levels, 2)
  for(i in 1:ncol(comparisons)){
    comp <- comparisons[,i]
    res[[i]] <-  tryCatch(pROC::roc(response = responseX,
                                        predictor = predictorX , levels=comp), error = function(x) NULL)
  }
  return(res)
}

#' Apply roc analysis on main factor on lowest level
#' @export
#' @examples
#' library(tidyverse)
#' library(SRMService)
#' config <- spectronautDIAData250_config$clone(deep=T)
#' config$parameter$min_nr_of_notNA  <- 20
#' data <- spectronautDIAData250_analysis
#' x <- sample(data$protein_Id,2)
#' data <- data %>% filter(protein_Id %in% x)
#' res <- compute_roc(data, config)
#' head(res)
#' i <- 1
#' pROC::plot.roc(res$rocs[[i]], print.auc = TRUE, main = paste(res$protein_Id[[i]], "\n",paste(res$rocs[[i]]$levels, collapse = " vs ")))
#' unique(res$protein_Id)
compute_roc <- function(data, config){
  nested <- data %>% group_by(!!sym(config$table$hierarchyKeys()[1]) ,
                              !!sym(config$table$hierarchyKeys(TRUE)[1])) %>% nest()
  nested <- nested %>% mutate(rocs = map(data ,
                                          rocs, response = config$table$factorKeys()[1],
                                          predictor= config$table$getWorkIntensity() ))

  nested <- nested %>% mutate(cls = map_lgl(rocs, is.null))  %>% dplyr::filter(cls == FALSE)
  #nested <- nested %>% mutate(names = map(rocs, names))

  dumm <- nested %>% dplyr::select(!!sym(config$table$hierarchyKeys()[1]),
                                   !!sym(config$table$hierarchyKeys(TRUE)[1]),
                                   rocs) %>%  unnest()
  head(dumm)
  dumm <- dumm %>% mutate(comparison = map_chr(rocs, function(x){paste(x$levels, collapse = " ")}))
  dumm <- dumm %>% separate(comparison, into = c("response1" , "response2"), sep=" ")
  dumm <- dumm %>% mutate(auc = map_dbl(rocs, pROC::auc)) %>% arrange(desc(auc))
  return(dumm)
}
