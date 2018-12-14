#
# Run script for WebGestaltR ORA analysis of MQ two group comparison on protein level
#

rm(list = ls())
library(limma)
library(SRMService)
library(WebGestaltR)
library(tidyverse)
library(quantable)

### Protein groups file
packagedir <- path.package("SRMService")

file = "yeast"

if(file == "yeast") {
  proteinGroupsFile <-
    file.path(packagedir, "/samples/proteinGroups/proteinGroupsYeast.txt")
  organism <- "scerevisiae"
  Experimentname <- "yeast_example"
  indx <- 4
} else if(file == "mouse") {
  proteinGroupsFile <-
    file.path(packagedir, "/samples/proteinGroups/proteinGroups.txt")
  organism <- "mmusculus"
  Experimentname <- "mouse_example"
  indx <- 3
} else if(file == "pullDown") {
  proteinGroupsFile <-
    file.path(packagedir, "/samples/proteinGroups/proteinGroupsPullDown.txt")
  organism <- "hsapiens"
  Experimentname <- "pullDown_example"
  indx <- 4
}


###


protein <- readr::read_tsv(proteinGroupsFile)
colnames(protein) <- make.names(colnames(protein))
tmp <- cumsum(rev(table(protein$Peptides)))
barplot(tmp[(length(tmp) - 5):length(tmp)], ylim = c(0, length(protein$Peptides)), xlab =
          'nr of proteins with at least # peptides')


rawF <-
  gsub("Intensity\\.", "", grep("Intensity\\.", colnames(protein), value =
                                  T))
condition <- quantable::split2table(rawF)[, indx]

annotation <- data.frame(
  Raw.file = rawF,
  Condition = condition,
  BioReplicate = paste0("X", 1:length(condition)),
  Run = 1:length(condition),
  IsotopeLabelType = rep("L", length(condition)),
  stringsAsFactors = F
)

###################################
### Configuration section
resultdir <- "output"
dir.create(resultdir)


nrNas = 5
nrPeptides = 2
reference = unique(annotation$Condition)[1]
qvalueThreshold = 0.05
qfoldchange = 1
numberOfProteinClusters = 3
enrichDatabase = "pathway_KEGG"
write.table(annotation, file = file.path(resultdir, "annotationused.txt"))


####### END of user configuration ##

# source("R/Grp2Analysis.R")
grp2 <- Grp2Analysis(
  annotation,
  Experimentname,
  maxNA = nrNas,
  nrPeptides = nrPeptides,
  reference = reference,
  numberOfProteinClusters = numberOfProteinClusters
)

grp2$setMQProteinGroups(protein)
grp2$setQValueThresholds(qvalue = qvalueThreshold, qfoldchange = qfoldchange)
mqQuantMatrixGRP2 <- grp2

head(mqQuantMatrixGRP2$getModPValuesCI())
usethis::use_data(mqQuantMatrixGRP2, overwrite = TRUE)

dir.create("output/ORA_inputFiles/")

tmp <- grp2$getNormalized()$data

ref_protein_list <- data.frame(IDs = row.names(tmp)) %>%
  separate(col = IDs,
           sep = "\\|",
           into = as.character(1:3)) %>%
  filter(`1` == "sp" | `1` == ">sp") %>%
  select(`2`)

write_tsv(ref_protein_list, "output/referencelist.txt", col_names = F)

clustering <-
  simpleheatmap3(t(tmp), labCol = row.names(tmp), plot = FALSE)$Col

clusterIDs <- clustering %>%
  group_by(clusterID) %>%
  separate(col = colID,
           sep = "\\|",
           into = as.character(1:3)) %>%
  filter(`1` == "sp" | `1` == ">sp") %>%
  select(clusterID, `2`) %>%
  ungroup()


write_interesting_geneFile <- function(ID, df, output.dir) {
  df %>%
    filter(clusterID == ID) %>%
    select(`2`) %>%
    write_tsv(
      path = paste(output.dir, "/protein_cluster_", ID, ".txt", sep = ""),
      col_names = F
    )
  return(paste("protein_cluster_",
               ID,
               ".txt was written to ",
               output.dir,
               sep = ""))
}

sapply(
  1:grp2$getNumberOfClusters(),
  write_interesting_geneFile,
  df = clusterIDs,
  output.dir = "output/ORA_inputFiles"
)

output <-
  WebGestaltR_batch(
    enrichMethod = "ORA",
    organism = organism,
    enrichDatabase = enrichDatabase,
    interestGeneFolder = "output/ORA_inputFiles",
    referenceGeneFile = "output/referencelist.txt",
    is.output = FALSE,
    interestGeneType = "uniprot_swissprot",
    referenceGeneType =  "uniprot_swissprot",
    outputDirectory = "output/"
  )

aggregate_list <- function(ll) {
  if (is.null(ll$enrichResult) || grepl("ERROR", ll$enrichResult))
    return(NULL)
  else {
    tmp <- as.data.frame(ll$enrichResult)
    tmp %>%
      mutate(file.origin = parse_number(tools::file_path_sans_ext(basename(ll$filename))))
  }
}

aggregated_results <- lapply(output, aggregate_list) %>%
  do.call(rbind, .)

# #REMOVE to render
rmarkdown::render(
  "vignettes/WebGestaltR_Grp2Analysis.Rmd",
  bookdown::html_document2(number_sections = FALSE),
  params = list(grp = mqQuantMatrixGRP2),
  clean = TRUE
)
