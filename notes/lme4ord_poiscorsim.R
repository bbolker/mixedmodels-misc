library(plyr)
library(lme4ord)

getCorPar <- function(x,unconstr=FALSE) {
    transEnv <- environment(x$parsedForm$random$NA.nlmeCorStruct$Lambdat$trans)
    corObj <- transEnv$object
    return(coef(corObj,unconstr=unconstr))
}    

simCor1 <- function(phi=0.8,sdgrp=2,sdres=1,
                    npergrp=20,ngrp=20,
                    seed=NULL,
                    ## set linkinv/simfun for GLMM sims
                    linkinv=identity,
                    simfun=identity) {
    if (!is.null(seed)) set.seed(seed)
    cmat <- sdres*phi^abs(outer(0:(npergrp-1),0:(npergrp-1),"-"))
    errs <- MASS::mvrnorm(ngrp,mu=rep(0,npergrp),Sigma=cmat)
    ranef <- rnorm(ngrp,mean=0,sd=sdgrp)
    d <- data.frame(f=rep(1:ngrp,each=npergrp))
    eta <- ranef[as.numeric(d$f)] + c(t(errs)) ## unpack errors by row
    mu <- linkinv(eta)
    d$y <- simfun(mu)
    d$tt <- factor(rep(1:npergrp,ngrp))
    return(d)
}


corObj <- NULL ## hack
corsimfun <- function(phi,sdgrp=2,sdres=1) {
    dp_sim <- simCor1(phi,sdgrp,sdres,
              linkinv=exp,simfun=function(x) rpois(length(x),lambda=x))
    corObj <<- nlme:::Initialize(nlme:::corAR1(0, form = ~ 1|f), dp_sim)
    r <- strucGlmer(y ~ 1 + (1|f)+nlmeCorStruct(1, corObj = corObj, sig = 1),
                    family=poisson,
                    data=dp_sim)
    res <- c(grpvar=c(VarCorr(r)[[1]]),ar1var=VarCorr(r)[[2]][1,1],
             phi=unname(getCorPar(r)))
    return(res)
}
set.seed(101)
## raply runs into trouble? some eval thing?
corvec <- seq(-0.9,0.9,by=0.1)
reslist <- setNames(vector(length(corvec),mode="list"),corvec)
for (i in seq_along(corvec)) {
    reslist[[i]] <- t(replicate(50,corsimfun(corvec[i])))
    saveRDS(reslist,file="lme4ord_poissim1.rds")
}

