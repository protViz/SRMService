#' extract intensities and annotations from MQ proteinGroups.txt
#' @export
#'
tidyMQ_ProteinGroups <- function(MQProteinGroups){
  pint <- select(MQProteinGroups, "Protein.group.Id" = "id", starts_with("Intensity."))
  pintLFQ <- select(MQProteinGroups, "Protein.group.Id" = "id", starts_with("LFQ.Intensity."))
  meta <- data.frame(ProteinName = MQProteinGroups$Majority.protein.IDs,
                     TopProteinName = sapply(strsplit(MQProteinGroups$Majority.protein.IDs, split=";"),
                                             function(x){x[1]}),
                     nrPeptides = MQProteinGroups$Peptides,
                     Fasta.headers = MQProteinGroups$Fasta.headers,
                     Protein.group.IDs = MQProteinGroups$id,
                     stringsAsFactors = F
  )
  pint <- pint %>%
    gather(key="Raw.file", value="MQ.Protein.Intensity", starts_with("Intensity.")) %>%
    mutate(Raw.file = gsub("Intensity.","",Raw.file))
  pintLFQ <- pintLFQ %>%
    gather(key="Raw.file", value="MQ.Protein.LFQ.intensity", starts_with("LFQ.intensity.")) %>%
    mutate(Raw.file = gsub("LFQ.intensity.","",Raw.file))

  pint <- inner_join(pint, pintLFQ , by=c("Protein.group.IDs","Raw.file"))
  res <- inner_join(meta, pint , by="Protein.group.IDs")
  return(res)
}


tidyMQ_Peptides <- function(MQPeptides){
  pint <- select(MQPeptides,"Peptides.Id"= "id", starts_with("Intensity."))
  idtype <- select(MQPeptides, "Peptides.Id"="id", starts_with("Identification.type."))
  meta <- select(MQPeptides, "Peptides.Id" = "id",
                 "Sequence",
                 "Proteins",
                 "Leading.razor.protein",
                 "Protein.group.IDs",
                 "Score",
                 "Missed.cleavages")

  PepIntensities <- pint %>%
    gather(key="Raw.file", value="Peptide.Intensity", starts_with("Intensity.")) %>%
    mutate(Raw.file = gsub("Intensity.","",Raw.file))


  PepIDType <- idtype %>%
    gather(key="Raw.file", value="IdType", starts_with("Identification.type.")) %>%
    mutate(Raw.file = gsub("Identification.type.","",Raw.file))

  tmp <-inner_join(PepIntensities,PepIDType, by=c("Peptides.Id", "Raw.file" ))
  xx<-inner_join(meta , tmp, by="Peptides.Id")

  xx$Proteotypic <-!grepl(";",xx$Protein.group.IDs)
  xx <- xx %>% separate_rows(Protein.group.IDs, sep=";",convert =TRUE)
  return(xx)
}



tidyMQ_Evidence <- function(Evidence){
  tmp <- select(Evidence,
                "Evidence.Id" = "id",
                "Peptide.Id" = "Peptide.ID",
                "Raw.file",
                "Protein.group.IDs",
                "Score",
                "Delta.score",
                "Calibrated.retention.time",
                "Charge",
                "MS.MS.count",
                "MS.MS.scan.number",
                "Evidence.Intensity" = "Intensity")
}
