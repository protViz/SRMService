rm(list=ls())
library(tidyverse)
library(readxl)
library(rlang)
library(yaml)
options(warn=0)

config <- createSpectronautPeptideConfiguration()
config$table$factors[["coding"]] = "coding"
config$table$factors[["sex"]] = "sex"
config$table$factors[["age"]] = "age"
config$table$factors[["Sample_id"]] = "Sample.Name"


#PAnnotated <- read.csv("SpectronautLibrary_Annotated_IncludingControls2.txt", sep="\t", stringsAsFactors = FALSE, dec=",")

PAnnotated <- SRMService::spectronautDIAData250
PAnnotated$Isotope.Label <- "light"
#head(PAnnotated %>% dplyr::filter(grepl("_LPQQANDYLNSFNWER_",PAnnotated$EG.ModifiedSequence)))
PAnnotated <- PAnnotated %>% dplyr::filter(!grepl("_Decoy$",PG.ProteinAccessions ))




resData <- setup_analysis(PAnnotated, config)


newcol <- paste("log2_", config$table$getWorkIntensity(), sep="")
resDataLog <- resData %>% mutate_at(config$table$workIntensity, .funs = funs(!!sym(newcol) := log2(.)))
config$table$setWorkIntensity(newcol)
config$parameter$workingIntensityTransform = "log"




# normalize data ----
wideMatrix <- toWideConfig(resDataLog, config, as.matrix = TRUE)
wideMatrix <- scale(wideMatrix)
resDataLog <- gatherItBack(wideMatrix,"srm_ScaledLogIntensity",config, resDataLog)

ggplot(resDataLog, aes_string(x = config$table$getWorkIntensity(), colour = config$table$sampleName )) +
  geom_line(stat="density")

proteins <- MSqRob::df2protdata(data.frame(resDataLog),
                                acc_col = "protein_Id",
                                run_name = "R.FileName",
                                quant_cols = "srm_ScaledLogIntensity",
                                annotations = NULL)

class(proteins)
slotNames(proteins)
slotNames(proteins@data[[1]])

#Fixed effects
fixed <- c("response")
#Random effects, for label-free data, it is best to always keep "Sequence"
random <- c("fragment_Id","R.FileName" )

#Do you want to save the model?
save_model <- TRUE #Alternative: save_model <- FALSE

#To which folder do you want to export the Excel file(s) and the saved model file (if you chose to do so)? Please do not use a trailing "/" here!
export_folder <- "./tmp/"

#Construct the contrast matrix L for testing on the fold changes of interest (i.e. our research hypotheses)
L <- makeContrast(contrasts=c("responseR - responseNR",
                              "responseR - responsePR",
                              "responsePR - responseNR",
                              "(responsePR + responseNR + responseR)/3-responsec"),
                  levels=c("responseR","responseNR","responsePR", "responsec"))

#Set the significance threshold (default: 5% FDR)
FDRlevel=0.05

#Fit the models
models <- fit.model(protdata=proteins, response="quant_value", fixed=fixed, random=random)

#Test the appropriate research hypotheses
results <- test.contrast_adjust(models, L)


names(results)[1]
for(i in 1:length(results)){
  tmp <- quantable::matrix_to_tibble(results[[i]])
  tmp$contrast <- names(results)[i]
  results[[i]] <- tmp
}

allContrasts <- bind_rows(results)
head(allContrasts)

write.csv(results,file=file.path(outdir,"MSqRobSignificances.tsv"))

head(allContrasts)
library(quantable)
pdf(file.path(outdir,"MSqRob_pValue_Volcano.pdf"))
quantable::multigroupVolcano(allContrasts,effect = "estimate", type="pval", condition="contrast",xintercept = c(-1,1),label = )
dev.off()

pdf(file.path(outdir,"MSqRob_qValue_Volcano.pdf"))
quantable::multigroupVolcano(allContrasts,effect = "estimate", type="qval", condition="contrast",xintercept = c(-1,1),label = )
dev.off()
