library(lme4)
library(blme)
library(glmmTMB)
library(brms)

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

riobglmer <- bglmer(form,
                    data=mydata, family=binomial (link='logit'),
                    fixef.prior= normal(cov = diag(9,12)))
rioglmer <- glmer(form,
                    data=mydata, family=binomial (link='logit'))
rioglmmTMB <- glmmTMB(form,
                    data=mydata, family=binomial (link='logit'))
prior <- get_prior(form,
                    data=mydata, family=bernoulli)
prior$prior[1] <- "normal(0,3)"
riobrms <- brm(contrast~c.con.tr*c.type.tr*c.diff.tr+(1|id)+(1|item.new),
               prior=prior,
               chains=3,
               data=mydata, family=bernoulli)
resList <- list(blme=riobglmer,lme4=rioglmer,glmmTMB=rioglmmTMB,
                brms=riobrms)
saveRDS(resList,file="bglmer_runs.rds")



