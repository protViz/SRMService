rm(list=ls())
library(reshape2)
library(plyr)
library(quantable)
longm <- read.csv("data/immunoLongFormat.txt",row.names = 1,stringsAsFactors = F)

longm$Intensity <- log2(longm$Intensity+1)
longm$Protein <- longm$Variable

data2Conditions <- longm

tmplong <- longm
nvariables <-3

data2Conditions <- longm
source("ProteinData.R")
protData <- Protein$new(data2Conditions)


data2Conditions$Condition[data2Conditions$Condition == "Healthy" | data2Conditions$Condition == "Gingivitis"] <- "NonPeriodontitis"
data2Conditions$Condition[(data2Conditions$Condition == "Aggressive" | data2Conditions$Condition == "Chronic")] <- "Periodontitis"


comparisonName <- "Periodontitis vs NonPeriodontitis"
rmarkdown::render("ROC2Group.Rmd", output_format = "html_document",
                  output_file = paste("ROC2GroupPeriodontitisVSNon.html",sep=""))
rmarkdown::render("ROC2Group.Rmd", output_format = "pdf_document",
                  output_file = paste("ROC2GroupPeriodontitisVSNon.pdf",sep=""))



Condition2Compare <- "Healthy|Gingivitis"
rmarkdown::render("SelectVars.Rmd",
                  output_file = "SelectVarsNonParandonditisVSParandontitis.html")
rmarkdown::render("SelectVars.Rmd", output_format = "pdf_document",
                  output_file = "SelectVarsNonParandonditisVSParandontitis.pdf")

