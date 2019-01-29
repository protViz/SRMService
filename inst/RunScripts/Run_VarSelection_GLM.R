rm(list=ls())
library(reshape2)
library(plyr)
library(dplyr)
library(quantable)
library(SRMService)



packagedir <- path.package("SRMService")
longm <- read.csv(file.path(packagedir, "samples/immunoLongFormat.txt"),row.names = 1,stringsAsFactors = F)

longm$Intensity <- log2(longm$Intensity+1)
longm$Protein <- longm$Variable


data2Conditions$Condition[data2Conditions$Condition == "Healthy" | data2Conditions$Condition == "Gingivitis"] <- "NonPeriodontitis"
data2Conditions$Condition[(data2Conditions$Condition == "Aggressive" | data2Conditions$Condition == "Chronic")] <- "Periodontitis"

unique(data2Conditions$Condition)

protData <- Protein$new(data2Conditions)


comparisonName <- "Periodontitis vs NonPeriodontitis"
rmarkdown::render("VariableSelection_ROCSingleProtein.Rmd", output_format = "html_document",
                  output_file = paste("ROC2GroupPeriodontitisVSNon.html",sep=""), clean = FALSE)
rmarkdown::render("VariableSelection_ROCSingleProtein.Rmd", output_format = "pdf_document",
                  output_file = paste("ROC2GroupPeriodontitisVSNon.pdf",sep="") , clean = FALSE)


tmplong <- longm
nvariables <-3


Condition2Compare <- "NonPeriodontitis"

rmarkdown::render("VariableSelection_GLM.Rmd",
                  output_file = "SelectVarsNonParandonditisVSParandontitis.html", clean = FALSE)

rmarkdown::render("VariableSelection_GLM.Rmd", output_format = "pdf_document",
                  output_file = "SelectVarsNonParandonditisVSParandontitis.pdf", clean=FALSE)

