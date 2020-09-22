library(tidyverse)
library(broom.mixed)
library(jglmm)
library(glmmTMB)
library(lme4)
library(GLMMadaptive)

jglmm_setup()

timefun <- function(x) {
    tm <- system.time(x <- eval(x,parent.frame()))
    attr(x,"time") <- tm
    return(x)
}

dd <- readRDS("gdata.rds")
ddx <- readRDS("gdata_disagg.rds")

savefun <- function() {
    nm <- ls(pattern="agg(_|$)", envir=parent.frame())
    print(nm)
    mod_list <- mget(nm, parent.frame())
    print(names(mod_list))
    save(mod_list, file="OH_batch.RData")
}

form_agg <- prop_succ~Intervention*Sex + (1|Subject) + (1|Image)
form_disagg <- update(form_agg, response ~ .)

glmer_agg <- timefun(glmer(form_agg, data= dd, family=binomial, weights=tot))
savefun()
glmer_disagg <- timefun(update(glmer_agg, form_disagg, data=ddx, weights=NULL))
savefun()
glmmTMB_agg <- timefun(glmmTMB(form_agg, data=dd, family=binomial, weights=tot))
savefun()
glmmTMB_disagg <- timefun(update(glmmTMB_agg, form_disagg, data=ddx, weights=NULL))
savefun()
glmmTMB_disagg_par2 <- timefun(update(glmmTMB_agg, form_disagg, data=ddx, weights=NULL,
                                      control=glmmTMBControl(parallel=2)))
savefun()
glmmTMB_disagg_nopar <- timefun(update(glmmTMB_agg, form_disagg, data=ddx, weights=NULL,
                                       control=glmmTMBControl(parallel=1)))
savefun()

## load("OH_batch.RData")
## McMasterPandemic::unpack(mod_list)

## dummy run for Julia warmup
print(timefun(jglmm(form_agg, data=dd, family="binomial", weights=dd$tot)))
jglmm_agg <- timefun(jglmm(form_agg, data=dd, family="binomial", weights=dd$tot))
savefun()
jglmm_disagg <- timefun(jglmm(form_disagg, data=ddx, family="binomial"))
savefun()


