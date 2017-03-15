rem shiny::runApp("inst/shiny/2Group2Test",port=1234, host="130.60.81.134")

"R.exe" -e "Sys.setenv(RSTUDIO_PANDOC='C:/Program Files/RStudio/bin/pandoc');shiny::runApp('C:/Users/wolski/prog/SRMService/inst/shiny/2Group2Test', port = 1234, host='130.60.81.134')"

rem will run 2Grp Shiny on server 130.60.81.134:1234



