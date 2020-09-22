## utilities/misc code for finding interesting packages related to mixed models

library(tidyverse)
library(miniCRAN)
library(crandep)
library(igraph)

rd <- c("Reverse depends", "Reverse imports", "Reverse linking to", "Reverse suggests")

dd <- get_dep_df("lme4",rd)
plot(igraph::graph_from_data_frame(dd))  ## too many!
## grep/regular expression tricks

a1 <- available.packages()
grep("lmm",rownames(a1),value=TRUE,ignore.case=TRUE)
## lm followed by ("m or e") followed by (a character not "t" or the end of the string)
grep("lm[me]([^t]|$)",rownames(a1),value=TRUE,ignore.case=TRUE)

