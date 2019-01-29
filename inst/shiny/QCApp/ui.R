## TESTESTTEST ##

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  # Application title
  titlePanel("QC Report"),
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      selectInput("selectFiletyp", label = h3("Select Filetype"),
                  choices = c("MaxQuant ProteinGroups" = "maxProt","MaxQuant Peptide.txt" = "maxPep",
                              "Progenesis Protein xls" = "progProt", "Progenesis Pepitde xls" = "progPep")
                  ,selected = 1),

      fileInput('proteinGroups', 'Choose MQ ProteinGroups File',
                accept=c('text/csv',
                         'text/comma-separated-values,text/plain',
                         '.csv','.zip')
                ),

      tags$hr(),
      uiOutput("parameterUI"),
      tags$hr(),
      uiOutput("generatereportbutton"),
      tags$hr(),
      uiOutput("downloadreport")
    ),
    # Show a plot of the generated distribution
    mainPanel(
      htmlOutput("fileInformation")
    )
  )
))
