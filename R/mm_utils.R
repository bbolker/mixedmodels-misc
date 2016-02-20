## miscellaneous mixed-model utilities

##' Wald confidence intervals (conditional on estimated theta
##' values)
easyPredCI <- function(model,newdata=NULL,alpha=0.05) {
    ## baseline prediction, on the linear predictor (logit) scale:
    pred0 <- predict(model,re.form=NA,newdata=newdata)
    ## fixed-effects model matrix for new data
    X <- model.matrix(formula(model,fixed.only=TRUE)[-2],newdata)
    beta <- fixef(model) ## fixed-effects coefficients
    V <- vcov(model)     ## variance-covariance matrix of beta
    pred.se <- sqrt(diag(X %*% V %*% t(X))) ## std errors of predictions
    ## inverse-link function
    linkinv <- family(model)$linkinv
    ## construct 95% Normal CIs on the link scale and
    ##  transform back to the response (probability) scale:
    crit <- -qnorm(alpha/2)
    linkinv(cbind(conf.low=pred0-crit*pred.se,
                  conf.high=pred0+crit*pred.se))
}

## Estimate 'denominator df' following the Pinheiro and Bates
## algorithm
calcDenDF <- function(fixed,grps,data,random=NULL,n0=1) {
    ## n0 is a hack (still not matching some values properly)
    ## try to construct random effects hierarchy
    ## see whether effects vary within groups
    testlevel1 <- function(grp,dmat) {
        vv <- apply(dmat,2,
                    function(x)
                    tapply(x,list(grp),function(x) length(unique(x))))
        apply(vv>1,2,any)
    }
    ##
    vlist <- list()
    if (!is.null(random)) {
        grps <- character(0)
        if (!is.list(random)) random <- list(random)
        lapply(random,
               function(x) {
            ## FIXME:: should do this via formula manipulation
                   v <- deparse(x[[2]][[2]]) ## left side of bar
                   g <- deparse(x[[2]][[3]]) ## grouping factor (NB may be wrong for complex groupings)
                   v2 <- setdiff(colnames(model.matrix(reformulate(v),data=data)),"(Intercept)")
                   vlist[[g]] <<- c(vlist[[g]],v2)
                   grps <<- union(grps,g)  ## NB may not be correctly ordered for multiple groups
               })
        ## hack to combine
        vnames <- unlist(Map(rep,as.list(names(vlist)),sapply(vlist,length)))
        vlist <- setNames(unlist(vlist),vnames)
    }
    ## get grouping factors even if they are expressed as interactions
    fdata <- lapply(data,factor)
    gfacs <- data.frame(setNames(lapply(grps,
                                        function(x) eval(parse(text=x),list2env(fdata))),grps),
                        check.names=FALSE)

    X <- model.matrix(fixed,data=data)
    np <- ncol(X)
    has.int <- "(Intercept)" %in% colnames(X)  ## FIXME: not tested
    ## apply level-variation function to each grouping factor
    lvary <- t(sapply(gfacs,testlevel1,dmat=X))
    ## add intercept and observation levels
    lvary <- rbind(pop=c(FALSE,rep(TRUE,np-1)),
                   lvary,
                   obs=rep(FALSE,np))
    ##
    for (i in seq_along(vlist)) {
        w <- which(names(vlist[i])==rownames(lvary))
        lvary[w,vlist[i]] <- FALSE
    }
    ## find index of first FALSE value 
    lev <- apply(lvary,2,function(x) which(!x)[1])
    p <- table(factor(lev,levels=seq(nrow(lvary))))  ## number of parameters estimated at each level
    N <- nrow(X)
    n <- c(pop=n0,sapply(gfacs,function(x) length(unique(x))),obs=N)  ## number of obs at each level
    ndiff <- n[-1] - n[-length(n)]
    dfs <- c(NA,setNames(ndiff-p[-1],names(p[-1])))
    r <- setNames(dfs[as.character(lev)],names(lev))
    if (has.int) {
        r["(Intercept)"] <- dfs[length(dfs)]
    } else {
        r <- r[-1]
    }
    r
}

## tests
## fixed <- ~poly(cyear,2)*sitetype
## grps <- "site"
## data <- avsiteX6
## calcDenDF(fixed,grps,data)

## calcDenDF(~age,"Subject",nlme::Orthodont)
## calcDenDF(~age,"Subject",nlme::Orthodont)
## calcDenDF(~age,data=nlme::Orthodont,random=~age|Subject)
## calcDenDF(~age,"Subject",nlme::Orthodont,random=~1|Subject)
## calcDenDF(~age,data=nlme::Orthodont,random=~1|Subject)

## intercept is wrong
## age ddf should be 26, not 25

## avsiteX6 <- read.csv("avsiteX6.csv")
## avsiteX6$cyear <- avsiteX6$year-min(avsiteX6$year)
## calcDenDF(~poly(cyear,2)*sitetype,"site",data=avsiteX6)
## calcDenDF(~poly(cyear,2)*sitetype,random=~poly(cyear,2)|site,data=avsiteX6)

KRmodcomp.drop1 <- function(model) {
    require("pbkrtest")
    m <- model.matrix(model)
    np <- ncol(m)-1
    orig_names <- colnames(m)
    old_data <- eval(getCall(model)$data) ## hope for the best
    krlist <- vector(ncol(m)-1,mode="list")
    reff <- sapply(findbars(formula(model)),deparse)
    for (i in seq(np+1)[-1]) {
        cat(i,orig_names[i],"\n")
        xx <- m[,c(-1,-i)] ## assume intercept
        if (np>1) {
            colnames(xx) <- paste0("b",seq(ncol(xx)))
            newvars <- colnames(xx)
        } else newvars <- "1"
        newform <- reformulate(c(newvars,reff),response=".")
        newdata <- cbind(xx,old_data)
        mod0 <- update(model,newform,data=newdata)
        krlist[[i-1]] <- KRmodcomp(model,mod0)
    }
    krlist
}

## see also: afex::mixed()

## options(warn=2)
## kk <- KRmodcomp.drop1(a1)
## fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
## kk2 <- KRmodcomp.drop1(fm1)
