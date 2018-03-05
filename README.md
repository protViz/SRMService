# SRMService

For more information about this project please visit the wproject wiki
<https://github.com/protviz/srmservice/wiki>


## Install

```{r}
library(devtools)
install_github(c( 'protViz/quantable', 'protViz/bibliospec))


pkgs <- c('gplots','missForest',  'qvalue','bookdown')
install.packages(pkgs)

source("https://bioconductor.org/biocLite.R")
biocLite("limma")

install_github(c('protViz/SRMService'))

```
