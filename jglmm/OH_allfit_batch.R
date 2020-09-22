library(tidyverse)
library(broom.mixed)
library(lme4)

dd <- readRDS("gdata.rds")
ddx <- readRDS("gdata_disagg.rds")
form_agg <- prop_succ~Intervention*Sex + (1|Subject) + (1|Image)
form_disagg <- update(form_agg, response ~ .)

load("OH_batch.RData")

aa <- allFit(mod_list[["glmer_disagg"]])

save("aa", file="OH_allfit_batch.rda")
