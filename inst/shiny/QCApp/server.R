## TESTESTTEST ##
## shiny::runApp("inst/shiny/2Group2Test",port=1234, host="130.60.81.134")
## shiny::runApp('C:/Users/wolski/prog/SRMService/inst/shiny/2Group2Test', port = 1234, host=)

library(shiny)
library(SRMService)
#library(rhandsontable)
options(shiny.maxRequestSize=30*1024^2)

mqPort<-function(input,v_upload_file){
  print(input$selectFiletyp)
  if(input$selectFiletyp == "maxProt"){
    protein <- read.csv(v_upload_file$filenam$datapath, sep="\t",stringsAsFactors = F, header=TRUE)
    v_upload_file$protein <- protein
    ## Prepare annotation table ###
    rawF <- gsub("Intensity\\.", "", grep("Intensity\\.",colnames(protein),value=T) )
    condition <- quantable::split2table(rawF)[,3]

    ## number of peptides plot ####
    nrPep <- cumsum(rev(table(protein$Peptides)))
    nrPeptidePlot <- renderPlot(barplot(nrPep[(length(nrPep)-5):length(nrPep)],
                                        ylim=c(0, length(protein$Peptides)),
                                        xlab='nr of proteins with at least # peptides'))
    v_upload_file$minPeptides <- max(protein$Peptides)

    ## number of NA's plot ###
    pint <- protein[,grep("Intensity\\.",colnames(protein))]
    pint[pint==0] <-NA

    pint2 <- pint[protein$Peptides >= 2,]
    nrNA <- apply(pint , 1, function(x){sum(is.na(x))})
    nrNA2 <- apply(pint2 , 1, function(x){sum(is.na(x))})

    naPlot <- renderPlot({
      par(mfrow=c(1,2))
      barplot((table(nrNA)),xlab="nr of NA's per protein (1 or more peptides)")
      barplot((table(nrNA2)),xlab="nr of NA's per protein (2 or more peptides)")
    })

    v_upload_file$maxNA <- ncol(pint)
    v_upload_file$maxMissing <- ncol(pint) - 4
    version <- help(package="SRMService")$info[[1]][4]
    ## prepare gui output
    list(nrPeptidePlot,
         naPlot,
         HTML(paste(v_upload_file$filenam[1], v_upload_file$filenam$datapath ,
                    "nr row: " ,
                    dim(protein)[1],
                    "nr columns: ",
                    dim(protein)[2],
                    "package Version : ",version
                    ,sep="<br/>")))
  }
}

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  grp2 <- NULL
  v_upload_file <- reactiveValues(data = NULL)
  v_download_links <- reactiveValues(filename= NULL)

  ### observes file upload
  filename <- observeEvent(input$proteinGroups,{
    v_upload_file$filenam <- input$proteinGroups
  })

  ## Create output button
  output$generatereportbutton <- renderUI({
    if(is.null(v_upload_file$filenam[1])){
      NULL
    }else{
      actionButton("generateReport", "Generate Report" )
    }
  })

  ##
  output$fileInformation <- renderUI({
    if(is.null(v_upload_file$filenam[1])){
      HTML("Please select: <br/> a) the type of the file. <br/> The file.")
    }else{
      mqPort(input,v_upload_file)
    }
  })

  output$parameterUI <- renderUI({
    if(is.null(v_upload_file$filenam[1])){
      NULL
    }else{
      list(
        numericInput("minPeptides", "Nr of Peptides per protein:", 2, max= v_upload_file$minPeptides),
        numericInput("maxMissing", "Maximum number of NAs: ",
                     value = v_upload_file$maxMissing,
                     min=0,
                     max = v_upload_file$maxNA),
        textInput("experimentID", "Experiment Title Name", "p2084_knockout"))
    }
  })

  ## Create some summary of the loaded data
  progress <- function(howmuch, detail){
    incProgress(howmuch, detail = detail)
  }

  ## react on GO button
  ## this method does all the computation
  generateReport <- eventReactive(input$generateReport, {
    #here will processing happen!
    if(is.null(v_upload_file$protein)){
      print("DUMM")
    }

    if(input$selectFiletyp == "maxProt"){}
    withProgress(message = 'Generating MaxQuant Protein Report', detail = "part 0", value = 0,
                 {
                   ### Rendering report
                   library(SRMService)
                   print(names(input))
                   annotation <- input$fileInformation
                   print("Annotation!")
                   print(annotation)


                   grp2 <- QCAnalysis(input$experimentID,
                                      maxNA = input$maxMissing,
                                      nrPeptides =input$minPeptides)

                   grp2$setMQProteinGroups(v_upload_file$protein)
                   incProgress(0.1, detail = paste("part", "Set up objects"))


                   tmpdir <- tempdir()
                   workdir <- file.path(tmpdir, gsub(" |:","_",date()))
                   rmdfile <- file.path( path.package("SRMService") , "/reports/Grp2Analysis.Rmd" )

                   if(!dir.create(workdir)){
                     stopApp(7)
                   }
                   rmdfile2run <- file.path(workdir ,"QCAnalysis.Rmd")
                   print(rmdfile2run)
                   if(!file.copy(rmdfile , rmdfile2run)){
                     stopApp(7)
                   }

                   rmarkdown::render(rmdfile2run,
                                     bookdown::pdf_document2())
                   incProgress(0.1, detail = paste("part", "Rendering"))
                   print(dir())
                   v_download_links$pdfReport <- file.path(workdir, "Grp2Analysis.pdf")

                   ### Writing p-values
                   write.table(grp2$getResultTable(), file=file.path(workdir,"pValues.csv"), quote=FALSE, sep = "\t", col.names=NA)
                   incProgress(0.1, detail = paste("part", "report"))
                   v_download_links$tsvTable <- file.path(workdir,"pValues.csv")
                 } ## end progress
    ) ## end with progress
    return(v_download_links$filename)
  }) ## end eventReactive


  output$downloadreport <- renderUI({
    files <- generateReport()

    downloads <- c(downloadReport="Download Report (.pdf)", downloadData = "Data (.xls)")
    ll <- list()

    for(i in 1:length(downloads)){
      ll[[i]]<-downloadButton(names(downloads)[i], label=downloads[[i]])
    }
    return(ll)
  })

  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$experimentID, "xls", sep = ".")
    },

    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      print(v_download_links$tsvTable)
      # Write to a file specified by the 'file' argument
      file.copy(v_download_links$tsvTable,file)
    }
  )
  output$downloadReport <- downloadHandler(
    filename = function() {
      paste(input$experimentID, "pdf", sep = ".")
    },

    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      print(v_download_links$pdfReport)
      # Write to a file specified by the 'file' argument
      file.copy(v_download_links$pdfReport, file)
    }
  )

})
