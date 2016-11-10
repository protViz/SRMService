#' fix scaffold file
#' @export
#'
fixScaffoldQuantReport <- function(filename){
  scftab <- read.csv(filename,sep="\t", stringsAsFactors = FALSE,skip=2)
  scftab <- scftab[-nrow(scftab),]
  tmp <- colnames(scftab)[2:ncol(scftab)]
  scftab <- scftab[,-ncol(scftab)]
  colnames(scftab) <- tmp

  ACC<-gsub(" \\(\\+[0-9]+\\)$","",scftab$Accession.Number)
  ACC<-gsub(" \\[[0-9]\\]$","",ACC)
  rownames(scftab) <- ACC

  return(scftab)
}
