## TESTESTTEST ##

library(shiny)
options(shiny.maxRequestSize=30*1024^2)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  v_upload_file <- reactiveValues(data = NULL)
  v_download_links <- reactiveValues(filename= NULL)

  filename <- observeEvent(input$proteinGroups,{
    print("hello")
    v_upload_file$filenam <- input$proteinGroups
  })

  output$button <- renderUI({
    if(is.null(v_upload_file$filenam[1])){
      tags$h4("welcome!")
    }else{
      actionButton("go", "GO" )
    }
  })

  randomVals <- eventReactive(input$go, {
    #here will processing happen!
    print("stupid")
    withProgress(message = 'Generating data', detail = "part 0", value = 0, {
      incProgress(0.3, detail = paste("part", "QC"))
      v_download_links$filename[[1]] <- "test.pdf"

      Sys.sleep(1)
      incProgress(0.3, detail = paste("part", "report"))
      v_download_links$filename[[2]] <- "dumbarton.pdf"

      Sys.sleep(2)
      incProgress(0.4, detail = paste("part", "all the crepp done"))

    })
  })

  output$downloadlink2 <- renderUI({
    if(length(v_download_links$filename)==0){
      tags$h4("welcome!")
    }else{
      res <- list()
      for(i in  1:length(v_download_links$filename)){
        res[[i]]<-tags$h4(v_download_links$filename[[i]])
      }
      return(res)
    }
  })

  output$downloadlinks <- renderUI({
    randomVals()
    tags$a(href = "url", "url")
    #list(tags$a(href = url, url),tags$br(), tags$a(href = url, url))
  })

})
