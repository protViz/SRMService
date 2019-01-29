allData <- read.csv( file="data/longFormat.txt",row.names = 1, stringsAsFactors = FALSE)
head(allData)

pool <- 1
#pool <- i
print(pool)
library(SRMService)

data <- allData[allData$pool==pool,]
head(data)
tmp <- SRMService::piwotPiw(data)
srms <- SRMService(data, paste("p1930Pool",pool,sep=""), qvalue=0.05)

srms$maxNAHeavy(30)
srms$maxNALight(40)
srms$maxNAFC(30)
rmarkdown::render("SRMPool2.Rmd",output_format = "pdf_document",
                  output_file = paste("SRMPool_",pool,".pdf",sep="")
)
rmarkdown::render("SRMPool2.Rmd",output_format = "html_document",
                  output_file = paste("SRMPool_",pool,".html",sep="")
)


rownames(protquant$data) <- split2table(rownames(protquant$data),split="\\|")[,3]

cc <- c("ZZ_FGCZCont0031", "A2MG_HUMAN", "ENOA_HUMAN", "HEMO_HUMAN", "FIBB_HUMAN")
cc <- gsub("ZZ_|_HUMAN","",cc)

prtable <- protquant
prtable$setHouseKeepers(cc)
rmarkdown::render("NormalizeByProteins.Rmd",output_format = "pdf_document",
                  output_file = paste("NormalizeByProteins_",pool,".pdf",sep=""))
