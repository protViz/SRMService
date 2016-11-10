fixScaffoldQuantReport <- function(filename){
  scftab <- read.csv(filename,sep="\t", stringsAsFactors = FALSE,skip=2)
  scftab <- scftab[-nrow(scftab),]
  tmp <- colnames(scftab)[2:ncol(scftab)]
  scftab <- scftab[,-ncol(scftab)]
  colnames(scftab) <- tmp
  return(scftab)
}
