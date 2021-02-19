# SRMService

For more information about this project please visit the wproject wiki
<https://github.com/protviz/srmservice/wiki>


## Install
For Debian 10
You will need to have installed:

```
sudo apt-get install libssl-dev pandoc pandoc-citeproc

# LaTeX
sudo apt-get install texlive-fonts-recommended texlive-latex-recommended texlive-latex-extra texlive-latex-base
```


```{r}
install.packages(c("bookdown", "conflicted", "corrplot", "dplyr", "forcats", "gridExtra", "heatmap3", "limma", "missForest", "pander", "plyr", "purrr", "quantable", "R6", "reshape2", "rlang", "rmarkdown", "scales", "shiny", "tibble", "tidyr", "tidyverse", "usethis"))
install.packages(c("tidyverse","usethis"))

install.packages("BiocManager")
BiocManager::install("limma")
```


```{r}
library(devtools)
install_github(c('protViz/SRMService'))

```


## Deployment on shiny-server




```{r}
library(devtools)
devtools::build_vignettes()
devtools::build()
```

Afterwards install with

```
R CMD INSTALL SRMService_*.tar.gz
```


