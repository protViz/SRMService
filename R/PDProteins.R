#' massage PD inputFiles.txt - remove columns and fix column names
#' @export
#' @param inputFiles read _InputFiles.txt file using the readr::read_tsv and pass it to inputFiles
#'
fix_PD_inputFiles <- function(inputFiles){
  colnames(inputFiles) <- make.names(colnames(inputFiles))
  inputFiles$Raw.File <- basename(inputFiles$File.Name)
  inputFiles <- inputFiles %>%
    filter(grepl("\\.raw", File.Name)) %>%
    mutate(Raw.File = gsub("\\.raw", "", Raw.File))
  inputFiles <- inputFiles %>% select(Study.File.ID, Raw.File)
  return(inputFiles)
}


#' remap pd proteins to format compatible with grp2Analysis reference class.
#' can be then set using set using the setProteins method
#' @export
#' @param proteins read _Proteins.txt file using readr::read_tsv and pass it to the proteins parameter
#' @param inputFiles read _InputFiles.txt file using the readr::read_tsv and pass it to inputFiles
#'
remap_PD_ProteinTable <- function(proteins, inputFiles){
  inputFiles <- fix_PD_inputFiles(inputFiles)
  colnames(proteins) <- make.names(colnames(proteins))
  proteins <- proteins %>%
    filter(Master == "IsMasterProtein") %>%
    filter(Exp.q.value.Combined < 0.01)


  proteinsX <- proteins %>% select(TopProteinName = Accession, nrPeptides = Number.of.Peptides,Fasta.headers=FASTA.Title.Lines, contains("Abundance.F") )
  colnames(proteinsX) <- gsub("\\.Sample","", colnames(proteinsX))
  colnames(proteinsX) <- gsub("Abundance\\.","Intensity.", colnames(proteinsX))

  pattern <- paste("\\.",inputFiles$Study.File.ID,"$",sep="")

  colnames(proteinsX) <- quantable::xxx_replace_xxx(colnames(proteinsX), pattern , paste(".", inputFiles$Raw.File , sep="") )
  pint <- proteinsX[, grepl("Intensity\\.", colnames(proteinsX))]
  proteinTable <- data.frame(ProteinName = proteinsX$TopProteinName,
                             TopProteinName = sapply(strsplit(proteinsX$TopProteinName,
                                                              split = ";"), function(x) {
                                                                x[1]
                                                              }),
                             nrPeptides = proteinsX$nrPeptides,
                             Fasta.headers = proteinsX$Fasta.headers,
                             pint, stringsAsFactors = F
  )
  return(proteinTable)
}
