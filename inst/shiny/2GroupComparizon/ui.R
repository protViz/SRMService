#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

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
      tags$hr(),
      fileInput('proteinGroups', 'Choose MQ ProteinGroups File',
                accept=c('text/csv',
                         'text/comma-separated-values,text/plain',
                         '.csv','.zip')),
      downloadButton('downloadReport', label="Download Report (.pdf)")
      ),
    # Show a plot of the generated distribution
    mainPanel(
       tableOutput('contents')
    )
  )
))
