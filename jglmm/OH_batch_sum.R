library(tidyverse)
library(broom.mixed)
## library(jglmm)
## jglmm_setup()
library(glmmTMB)
load("OH_batch.RData")

## drop formulas and Julia models (serializing/deserializing is ugly)
mod_list <- mod_list[!grepl("^form|jglmm",names(mod_list))]
## tidying big glmmTMB models is slooooow ... ?
## tidy(mod_list[[4]], effects="fixed")
mod_sum <- map_dfr(mod_list,tidy, effects="fixed", .id="pkg", conf.int=TRUE)
mod_times <- map_dfr(mod_list,
                     ~ tibble(time=attr(.,"time")[["elapsed"]]),
                     .id="pkg")

jterms_to_Rterms <- function(tt)  {
    tt %>% str_remove_all("(:| )") %>% str_replace("&",":")
}

## get julia csv
tmp_tidy <- function(fn="julia_agg.csv", conf.level=0.95) {
    qq <-  qnorm((1+conf.level)/2)
    (read_csv(fn)
        %>% rename(estimate=`Estimate`,
                   std.error=`Std.Error`,
                   statistic=`z value`,
                   p.value=`P(>|z|)`)
        %>% mutate_at("term", jterms_to_Rterms)
        %>% mutate(conf.low=estimate-qq*std.error,
                   conf.high=estimate+qq*std.error)
    )
}

m2 <- map_dfr(list(jglmm_agg="julia_agg.csv", jglmm_disagg="julia_disagg.csv"),
              tmp_tidy, .id="pkg")
oh_mod_sum <- bind_rows(mod_sum,m2) %>% select(-c(effect,statistic,p.value,component))

jt <- tibble(pkg=c("jglmm_agg","jglmm_disagg"),time=c(0.688,77.510714))
oh_mod_times <- bind_rows(mod_times,jt)

save(oh_mod_sum, oh_mod_times, file="OH_batch_sum.rda")
