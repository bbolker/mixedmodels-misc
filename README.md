## mixed models miscellany

A variety of mixed-model-related resources.

### Task view

A (so far very preliminary) draft of a [Task View](http://cran.r-project.org/web/views/) for mixed models in R.

The format is described in the [Writing CRAN Task Views vignette](http://cran.r-project.org/web/packages/ctv/vignettes/ctv-howto.pdf) of the [ctv package](http://cran.r-project.org/web/packages/ctv/index.html).

To view the task view as HTML, from within R:
```
library("ctv")
mm_ctv <- read.ctv("MixedModels.ctv")
ctv2html(mm_ctv)
browseURL("./MixedModels.html")
```

### FAQ

`glmmFAQ`: An updated version of the generalized linear mixed models FAQ currently hosted at [wikidot.com](http://glmm.wikidot.com/faq) ([HTML](https://htmlpreview.github.io/?https://raw.githubusercontent.com/bbolker/mixedmodels-misc/master/glmmFAQ.html), [source](glmmFAQ.rmd))


### GLMM worked examples

`ecostats_chap`: examples from chapter 13  *Ecological Statistics: Contemporary Theory and Application*  editors Negrete, Sosa, and Fox, plus data ([HTML](https://htmlpreview.github.io/?https://raw.githubusercontent.com/bbolker/mixedmodels-misc/master/ecostats_chap.html), [source](ecostats_chap.rmd), [data](data))

### Various notes

- [thoughts about interfaces for mixed modeling software](mixed_interface.rmd)
- [thoughts about design of DSLs for statistics](stats_design.md)
