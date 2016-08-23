## TESTESTTEST ##

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

  # Application title
  titlePanel("2 Group Comparizon"),
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      numericInput("peptides", "Nr of Peptides per protein:", 2),
      numericInput("maxMissing", "Maximum number of missing values",8),
      textInput("experimentID", "Experiment ID","p2084"),
      tags$hr(),
      fileInput('proteinGroups', 'Choose MQ ProteinGroups File',
                accept=c('text/csv',
                         'text/comma-separated-values,text/plain',
                         '.csv','.zip')),
      uiOutput("button")


    ),
    # Show a plot of the generated distribution
    mainPanel(
      uiOutput("downloadlinks"),
      tags$hr(),
      uiOutput("downloadlink2")
    )
  )
))
