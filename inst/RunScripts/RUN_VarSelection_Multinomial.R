rm(list=ls())

library(SRMService)
library(reshape2)

packagedir <- path.package("SRMService")
longm <- read.csv(file.path(packagedir, "samples/immunoLongFormat.txt"),row.names = 1,stringsAsFactors = F)
longm$Intensity <- log2(longm$Intensity+1)



#head(datam)
longm$Protein <- longm$Variable
protData <- ProteinVariableSelect$new(longm)

rmarkdown::render("VariableSelection_MultinomialLasso.Rmd",output_format = "pdf_document", clean=FALSE)



