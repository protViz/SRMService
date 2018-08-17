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


## Deployment on server

To install vignettes with the package you first need to build the package with package vignettes. Only then you can install the package.

Building the package on the system you want to deploy it might fail. Therefore, build it first on the development system.

```{r}
library(devtools)
devtools::build()
```

On github create release and attach product.

On server get the product with `wget` and `R CMD INSTALL SRMService_*.tar.gz`
