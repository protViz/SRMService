#' extract intensities and annotations from MQ proteinGroups.txt
#' @export
MQProteinGroups <- function(MQProteinGroups){
  pint <- MQProteinGroups[,grep("Intensity\\.",colnames(MQProteinGroups))]
  proteinTable <- data.frame(ProteinName = MQProteinGroups$Majority.protein.IDs,
                             TopProteinName = sapply(strsplit(MQProteinGroups$Majority.protein.IDs, split=";"),
                                                     function(x){x[1]}),
                             nrPeptides = MQProteinGroups$Peptides,
                             Fasta.headers = MQProteinGroups$Fasta.headers,
                             pint, stringsAsFactors = F)
}
