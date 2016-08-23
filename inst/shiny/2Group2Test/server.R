## TESTESTTEST ##

library(shiny)
options(shiny.maxRequestSize=30*1024^2)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  v_upload_file <- reactiveValues(data = NULL)
  v_download_links <- reactiveValues(filename= NULL)


  ### observes file upload
  filename <- observeEvent(input$proteinGroups,{
    print("hello")
    v_upload_file$filenam <- input$proteinGroups
  })

  ## Create output button
  output$button <- renderUI({
    if(is.null(v_upload_file$filenam[1])){
      NULL
    }else{
      actionButton("go", "Generate Report" )
    }
  })

  ##
  output$fileInformation <- renderUI({
    if(is.null(v_upload_file$filenam[1])){
      ("Please upload a ProteinGroups.txt file.")
    }else{
      HTML(paste("Some brief summary of the file goes here.",
      "Therea 12 samples and 3200 proteins with at least 2 peptides",
      "max intensity, min intensity, nr of missing values per row",
      "maybee some plot to help to set the right number of missing values per protein.",
      "Than you can trigger the report generation with the Generate Report button."
      ,sep="<br/>"))

    }
  })
  ## Create some summary of the loaded data


  ## react on GO button
  ## this method does all the computation
  randomVals <- eventReactive(input$go, {
    #here will processing happen!

    withProgress(message = 'Generating data', detail = "part 0", value = 0, {
      incProgress(0.3, detail = paste("part", "QC"))
      v_download_links$filename[[1]] <- "test.pdf"

      Sys.sleep(1)
      incProgress(0.3, detail = paste("part", "report"))
      v_download_links$filename[[2]] <- "dumbarton.pdf"

      Sys.sleep(2)
      incProgress(0.4, detail = paste("part", "all the crepp done"))

    })
    return("will return list of file locations")
  })

  # output$downloadlinks <- renderUI({
  #   if(length(v_download_links$filename)==0){
  #     NULL
  #   }else{
  #     res <- list()
  #     for(i in  1:length(v_download_links$filename)){
  #       res[[i]]<-tags$h4(v_download_links$filename[[i]])
  #     }
  #     return(res)
  #   }
  # })

  output$downloadlink2 <- renderUI({
    files <- randomVals()

    downloads <- c(downloadReport="Download Report (.pdf)", downloadData = "Data (.tsv)")
    ll <- list()

    for(i in 1:length(downloads)){
      ll[[i]]<-downloadButton(names(downloads)[i], label=downloads[[i]])
    }
    return(ll)
  })

})
