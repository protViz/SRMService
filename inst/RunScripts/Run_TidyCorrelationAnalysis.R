rm(list=ls())
knitr::opts_chunk$set(echo = TRUE, message = FALSE, fig.width = 8, fig.height = 12)

library(conflicted)
library(tidyverse)
library(rlang)
library(dplyr)
library(ggplot2)

rm(list=ls())
library(SRMService)
createSpectronautPeptideConfiguration <- function(isotopeLabel="Isotope.Label", qValue="EG.Qvalue"){
  atable <- AnalysisTableAnnotation$new()
  atable$fileName = "R.FileName"

  # measurement levels.
  atable$hierarchy[["protein_Id"]]    <-  "PG.ProteinAccessions"
  atable$hierarchy[["peptide_Id"]]    <-  "PEP.StrippedSequence"
  atable$hierarchy[["modPeptide_Id"]] <-  "EG.ModifiedSequence"
  atable$hierarchy[["fragment_Id"]]   <-  c("EG.ModifiedSequence", "FG.Charge")

  #
  atable$qValue = qValue
  atable$startIntensity = "FG.Quantity"
  atable$workIntensity = "FG.Quantity"
  atable$isotopeLabel = isotopeLabel
  anaparam <- AnalysisParameters$new()
  configuration <- AnalysisConfiguration$new(atable, anaparam)
}

config <- createSpectronautPeptideConfiguration()
config$table$factors[["coding"]] = "coding"
config$table$factors[["sex"]] = "sex"
config$table$factors[["age"]] = "age"
config$table$factors[["Sample_id"]] = "Sample.Name"

config$table$hierarchy
R6extractValues(config)

if(0){
  not <- c("R.Condition", "R.Fraction", "R.Label", "R.Replicate", "PG.ProteinGroups" ,"PG.Cscore","PG.Qvalue","PG.RunEvidenceCount",
           "PEP.GroupingKey", "PEP.GroupingKeyType",   "PEP.IsProteotypic",               "PEP.NrOfMissedCleavages",
           "PEP.Rank",   "PEP.RunEvidenceCount",            "PEP.UsedForProteinGroupQuantity",
           "EG.iRTPredicted",                 "EG.IsDecoy" ,"EG.ModifiedPeptide",
           "EG.UserGroup",                    "EG.Workflow",
           "EG.IsUserPeak",                   "EG.IsVerified",  "EG.iRTEmpirical",
           "EG.MeanApexRT",                   "EG.RTPredicted",                  "EG.AvgProfileQvalue",             "EG.MaxProfileQvalue",
           "EG.MinProfileQvalue",             "EG.PercentileQvalue",             "EG.ReferenceQuantity..Settings.", "EG.TargetQuantity..Settings.",
           "FG.PrecMz",                       "FG.PrecMzCalibrated",             "FG.ShapeQualityScore",
           "Tube.Id","Extract.Name" ,"Species" ,"Condition","Resource_id", "EG.good")



  longFormat <- readr::read_tsv("D:/projects/p2244_ProjectOrder_3687/inputData/preprocessing/withAnnotation_3687_DB_WithIsoforms.txt")
  colnames(longFormat) <- make.names(colnames(longFormat))
  longFormat <- longFormat %>% select(which(!colnames(longFormat) %in% not))
  spectronautDIAData250 <- longFormat
  usethis::use_data(spectronautDIAData250,overwrite = TRUE)
}


longFormat <- get(data(spectronautDIAData250))
"PEP.StrippedSequence" %in% colnames(longFormat)
longFormat$Isotope.Label <- "Light"

#source("c:/Users/wolski/prog/SRMService/inst/RPrototypes/tidyMS_R6_TransitionCorrelations.R")

config$table$workIntensity <- config$table$workIntensity[1]
longNoDecoy <- setup_analysis(longFormat, config)


data <- longNoDecoy
longNoDecoy <- setIntensitiesToNA(longNoDecoy, config)

## QValue Summaries

longQSummaries <- summariseQValues(longNoDecoy, config)

SRMService::getMissingStats(longQSummaries, config)
longQNASummaries <- summariseNAs(longQSummaries, config)



## Filter

qvalFilt <- longQNASummaries %>% filter_at( "srm_QValueMin" , all_vars(. < config$parameter$qValThreshold )   )
SRMService::hierarchyCounts(qvalFilt, config)


tmp <- toWideConfig(qvalFilt, config)
tmp2 <- as.matrix(tmp[(length(config$table$hierarchyKeys() ) + 1):ncol(tmp)])
rownames(tmp2) <- tmp %>% unite(x ,1,4) %>% pull(x)

#quantable::simpleheatmap(cor(log2(tmp2), use="pairwise.complete.obs"), margin=c(15,5))


## Remove single hit wonders
# Drop all proteins with only one precursor.


#```{r}
qvalFilt2 <- nr_B_in_A(qvalFilt,config$table$hierarchyKeys()[1], config$table$hierarchyKeys()[2])
table(qvalFilt2$nr_peptide_Id_by_protein_Id)
qvalFiltV <- dplyr::filter(qvalFilt2, nr_peptide_Id_by_protein_Id > 1)
hierarchyCounts(qvalFiltV,config)

#```

# Correlation Filtering

#```{r}
options(warn=0)
xx <- markDecorrelated(qvalFiltV, config)
hierarchyCounts(xx, config)


qvalFiltCorr <- dplyr::filter(xx, srm_decorelated == FALSE)
hierarchyCounts(qvalFiltCorr, config)
intensitiesD <- toWideConfig(qvalFiltCorr, config)

#```



#```{r}
#quantable::simpleheatmap(cor(log2( intensitiesD[,sapply(intensitiesD, class) == "numeric"]),
#                             use="pairwise.complete.obs"), margin=c(15,5))

#```

### Rank precursors by NA

qvalFiltCorrS <- rankPrecursorsByNAs(qvalFiltCorr, config)
qvalFiltCorrS %>% select(protein_Id,fragment_Id, srm_NrNAs, srm_NrNARank) %>% distinct()
qvalFiltCorrSI <- rankPrecursorsByIntensity(qvalFiltCorrS, config)
config$table$getWorkIntensity()
config$table$workIntensity <- config$table$workIntensity[-length(config$table$workIntensity)]
qvalFiltCorrSI %>% select(protein_Id,fragment_Id, starts_with("srm_")) %>% distinct()

qvalFiltImputed <- impute_correlationBased(qvalFiltCorrSI, config)
qvalFiltImputed %>% select(protein_Id,fragment_Id, starts_with("srm_")) %>% distinct()

# Aggregation

# Rank precursors by intensity in all samples

#```{r}
#qvalFiltImputed <- rankPrecursorsByIntensity(qvalFiltImputed, config)

#```



#```{r}
proteinIntensities <- aggregateTopNIntensities(qvalFiltImputed,config,N=3)
#qvalFiltImputed %>% select(config$table$hierarchyKeys()[1], nr_peptide_Id_by_protein_Id) %>% arrange(nr_peptide_Id_by_protein_Id)

proteinsWide <- spread(proteinIntensities, sampleName,srm_sumTopInt)

#quantable::simpleheatmap(t(scale(t(log2(1+proteinsWide[2:ncol(proteinsWide)])))), margins = c(16,5))
nrPep <-qvalFiltImputed %>% select(config$table$hierarchyKeys()[1], nr_peptide_Id_by_protein_Id) %>% distinct()
xx <- inner_join(nrPep,proteinsWide)
#readr::write_csv(xx, path="output/proteinQuantsTidy.txt")

#```



