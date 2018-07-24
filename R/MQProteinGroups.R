#' extract intensities and annotations from MQ proteinGroups.txt
#' @export
#' @param MQProteinGroups data.frame generated with read.csv("peptide.txt",sep="\\t", stringsAsFactors=FALSE)
#' @examples
#' protein_txt <- system.file("samples/maxquant_txt/MSQC1/proteinGroups.txt",package = "SRMService")
#' protein_txt
#' protein_txt <- read.csv(protein_txt, header=TRUE, stringsAsFactors = FALSE, sep="\t")
#' res <-tidyMQ_ProteinGroups(protein_txt)
#' head(res)
tidyMQ_ProteinGroups <- function(MQProteinGroups){
  if(is.character(MQProteinGroups)){
    MQProteinGroups <- read.csv(MQProteinGroups, header=TRUE, stringsAsFactors = FALSE, sep="\t")
  }
  colnames(MQProteinGroups) <- tolower(colnames(MQProteinGroups))

  pint <- select(MQProteinGroups, "protein.group.id" = "id", starts_with("intensity."))
  pintLFQ <- select(MQProteinGroups, "protein.group.id" = "id", starts_with("lfq.intensity."))
  meta <- data.frame(protein.names = MQProteinGroups$majority.protein.ids,
                     top.protein.name = sapply(strsplit(MQProteinGroups$majority.protein.ids, split=";"),
                                             function(x){x[1]}),
                     nr.peptides = MQProteinGroups$peptides,
                     fasta.headers = MQProteinGroups$fasta.headers,
                     protein.group.id = MQProteinGroups$id,
                     stringsAsFactors = F
  )
  pint <- pint %>%
    gather(key="raw.file", value="mq.protein.intensity", starts_with("intensity.")) %>%
    mutate(raw.file = gsub("intensity.","",raw.file))

  pintLFQ <- pintLFQ %>%
    gather(key="raw.file", value="mq.protein.lfq.intensity", starts_with("lfq.intensity.")) %>%
    mutate(raw.file = gsub("lfq.intensity.","",raw.file))

  pint <- inner_join(pint, pintLFQ , by=c("protein.group.id","raw.file"))
  res <- inner_join(meta, pint , by="protein.group.id")
  return(res)
}

#' parse MQ peptides.txt
#' @export
#' @param MQPeptides data.frame generated with read.csv("peptide.txt",sep="\\t", stringsAsFactors=FALSE)
#' @examples
#' library(tidyverse)
#' peptides_txt <- system.file("samples/maxquant_txt/MSQC1/peptides.txt",package = "SRMService")
#' peptides_txt <- read.csv(peptides_txt, header=TRUE, stringsAsFactors = FALSE, sep="\t")
#' str(peptides_txt)
#' tmp <-paste(peptides_txt$Evidence.IDs, collapse = ";")
#' tmp <- strsplit(tmp, ";")
#' length(unique(tmp[[1]]))
#' res <-tidyMQ_Peptides(peptides_txt)
#' dim(res)
#' head(res)
tidyMQ_Peptides <- function(MQPeptides){
  if(is.character(MQPeptides)){
    MQPeptides <- read.csv(MQPeptides, header=TRUE, stringsAsFactors = FALSE, sep="\t")
  }
  colnames(MQPeptides) <- tolower(colnames(MQPeptides))
  pint <- select(MQPeptides,"peptides.id"= "id", starts_with("intensity."))
  idtype <- select(MQPeptides, "peptides.id"="id", starts_with("identification.type."))
  meta <- select(MQPeptides, "peptides.id" = "id",
                 "sequence",
                 "proteins",
                 "leading.razor.protein",
                 "protein.group.ids",
                 "score",
                 "pep",
                 "missed.cleavages")

  PepIntensities <- pint %>%
    gather(key="raw.file", value="peptide.intensity", starts_with("intensity.")) %>%
    mutate(raw.file = gsub("intensity.","",raw.file))


  PepIDType <- idtype %>%
    gather(key="raw.file", value="id.type", starts_with("identification.type.")) %>%
    mutate(raw.file = gsub("identification.type.","",raw.file))

  tmp <-inner_join(PepIntensities,PepIDType, by=c("peptides.id", "raw.file" ))
  xx<-inner_join(meta , tmp, by="peptides.id")

  xx$proteotypic <-!grepl(";",xx$protein.group.ids)
  xx <- xx %>% separate_rows(protein.group.ids, sep=";",convert =TRUE)
  return(xx)
}
#' read evidence file
#' @export
#' @examples
#' library(tidyverse)
#' evidence_txt <- system.file("samples/maxquant_txt/MSQC1/evidence.txt",package = "SRMService")
#' evidence_txt <- read.csv(evidence_txt, header=TRUE, stringsAsFactors = FALSE, sep="\t")
#' xx <- tidyMQ_Evidence(evidence_txt)
tidyMQ_Evidence <- function(Evidence){
  colnames(Evidence) <- tolower(colnames(Evidence))
  tmp <- select(Evidence,
                "evidence.id" = "id",
                "peptide.id",
                "raw.file",
                "protein.group.ids",
                "score",
                "delta.score",
                "calibrated.retention.time",
                "charge",
                "mass",
                "ms.ms.count",
                "ms.ms.scan.number",
                "evidence.intensity" = "intensity")
  return(tmp)
}
