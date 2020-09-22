library(glmmTMB)
library(lme4)
library(tidyverse)
library(drake)
library(broom.mixed)
library(jglmm)
jglmm_setup()
## FIXME: for timing purposes,
## run trivial model to JIT mixed model machinery?



## timing function (might be able to magically get this from drake info, but ???)
timefun <- function(x) {
    tm <- system.time(x <- eval(x,parent.frame()))
    attr(x,"time") <- tm
    return(x)
}

gfit_plan <- drake_plan(
    
    dd=(read_csv(file_in("OH_data.csv")) %>%
       mutate(tot=Successes+Failures,
              prop_succ=Successes/tot,
              Intervention=relevel(factor(Intervention),"Pre"),
              Sex=factor(Sex))
    ),
    
    ddx=(dd
        %>% group_by(Subject,Image,Intervention,Sex)
        %>% summarise(response=rep(c(1,0),c(Successes,Failures)),.groups="drop")
    ),
    
    form_agg=prop_succ~Intervention*Sex + (1|Subject) + (1|Image),
    
    form_disagg=update(form_agg, response ~ .),
    
    ## readd(): see https://github.com/ropensci/drake/issues/1163
    glmer_agg=timefun(glmer(form_agg, data=readd(dd), family=binomial, weights=tot)),
    
    glmer_disagg=timefun(update(glmer_agg, form_disagg, data=readd(ddx), weights=NULL)),
    
    glmmTMB_agg=timefun(glmmTMB(form_agg, data=readd(dd), family=binomial, weights=tot)),
    
    glmmTMB_disagg=timefun(update(glmmTMB_agg, form_disagg, data=readd(ddx), weights=NULL)),

    jglmm_agg=timefun(jglmm(form_agg, data=readd(dd), family="binomial", weights=readd(dd)$tot)) ,

    ## jglmm_disagg=timefun(jglmm(form_disagg, data=readd(ddx), family="binomial")) ,

    ## FIXME: is there a way to the equivalent of ls(pattern="agg$")? or deal with
    ##  static branching, transforms, etc. ???
    modList={
        lme4:::namedList(
                   glmer_agg,glmer_disagg, glmmTMB_agg, glmmTMB_disagg, jglmm_agg
                   ## , jglmm_disagg  ## skip for now
               )
    },
    
    modsum=(map_dfr(modList, ~tidy(., effects="fixed", conf.int=TRUE),.id="model")
            %>% separate(model,into=c("pkg","type"),sep="_",remove=FALSE)),
    
    out=saveRDS(modList,file_out("gfits.rds")),
    
    out_sum=saveRDS(modsum,file_out("gsum.rds")),
    
    out_data=saveRDS(dd,file_out("gdata.rds")),
    out_data_disagg=saveRDS(ddx,file_out("gdata_disagg.rds")),
    out_csv=write.csv(dd,file_out("gdata.csv")),
    out_disagg_csv=write.csv(ddx,file_out("gdata_disagg.csv"))
)
drake_config(gfit_plan)
