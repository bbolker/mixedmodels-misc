library(lme4)
library(blme)
library(glmmTMB)
library(brms)
library(MCMCglmm)

## @knitr setup_runs
form <- contrast~c.con.tr*c.type.tr*c.diff.tr+(1|id)+(1|item.new)
mydata <- expand.grid(c.con.tr=factor(1:3),
                      c.type.tr=factor(1:2),
                      c.diff.tr=factor(1:2),
                      id=factor(1:10),
                      item.new=factor(1:10))

## @knitr simulate_data
set.seed(101)
mydata$contrast <- simulate(form[-2],
                            newdata=mydata,
                            newparams=list(beta=rep(1,12),
                                           theta=rep(1,2)),
                            family=binomial,
                            weights=rep(1,nrow(mydata)))[[1]]

## @knitr run_models

t.bglmer <- system.time(fit.bglmer <- bglmer(form,
                    data=mydata, family=binomial (link='logit'),
                    fixef.prior= normal(cov = diag(9,12))))
t.glmer <- system.time(fit.glmer <- glmer(form,
                    data=mydata, family=binomial (link='logit')))
t.glmmTMB <- system.time(fit.glmmTMB <- glmmTMB(form,
                    data=mydata, family=binomial (link='logit')))
prior <- get_prior(form,
                    data=mydata, family=bernoulli)
prior$prior[1] <- "normal(0,3)"
t.brms <- system.time(fit.brms <- brm(contrast~c.con.tr*c.type.tr*c.diff.tr+(1|id)+(1|item.new),
               prior=prior,
               chains=3,
               data=mydata, family=bernoulli))

t.MCMCglmm <- system.time(fit.MCMCglmm <- MCMCglmm(contrast~c.con.tr*c.type.tr*c.diff.tr,
                        random=~id+item.new,
                        data=mydata,
                        prior=list(B=list(mu=rep(0,12),V=diag(9,12)),
                                   R=list(V=1,nu=0),
                                   G=list(list(V=1,nu=0),
                                          list(V1=1,nu=0)))))
resList <- list(blme=fit.bglmer,lme4=fit.glmer,glmmTMB=fit.glmmTMB,
                brms=fit.brms,MCMCglmm=fit.MCMCglmm)
attr(resList,"times") <-
    list(t.bglmer,t.glmer,t.glmmTMB,t.brms,t.MCMCglmm)
saveRDS(resList,file="bglmer_runs.rds")



