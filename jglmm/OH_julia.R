library(tidyverse)
library(jglmm)
library(broom.mixed)
dd <- readRDS("gdata.rds")

jglmm_setup()
form_agg <- prop_succ~Intervention*Sex + (1|Subject) + (1|Image)
jglmm_agg <- jglmm(form_agg, data=dd, family="binomial", weights=dd$tot)
jglmm(form_agg, data=dd, family="binomial", weights=dd$tot,
                   return_val="julia_model_str")
tidy(jglmm_agg, effects="fixed", conf.int=TRUE)
