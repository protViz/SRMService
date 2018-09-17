# SRMService

For more information about this project please visit the wproject wiki
<https://github.com/protviz/srmservice/wiki>


## Install

```{r}
library(devtools)
install_github(c('protViz/SRMService'))

```


## Deployment on server


```{r}
library(devtools)
devtools::build_vignettes()
devtools::build()
```

Afterwards install with

```
R CMD INSTALL SRMService_*.tar.gz
```


