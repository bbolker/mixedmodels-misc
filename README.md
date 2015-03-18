## mixed models task view

This is a (so far very preliminary) draft of a [Task View](http://cran.r-project.org/web/views/) for mixed models in R.

The format is described in the [Writing CRAN Task Views vignette](http://cran.r-project.org/web/packages/ctv/vignettes/ctv-howto.pdf) of the [ctv package](http://cran.r-project.org/web/packages/ctv/index.html).

To view the task view as HTML, from within R:
```
library("ctv")
x <- read.ctv("MixedModels.ctv")
ctv2html(x)
browseURL("./MixedModels.html")
```
