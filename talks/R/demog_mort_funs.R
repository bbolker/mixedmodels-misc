plotfun <- function(model,dat,
                    formula=~PERIOD*Previous.size*Treatment+Reef.size,
                    drop.cols=numeric(0),
                    colors=c("purple","gray"),
                    lwd=1,jitter.y=0.025,exp.y=0.08,title="",
                    incl.var=FALSE,
                    ctrval=1, ## centering value for sizes
                    ctrline=TRUE ## draw line at center value?
                    ) {
  ## prediction data frame
  mortframe <- with(dat,expand.grid(Treatment=levels(Treatment),
                                    PERIOD=levels(PERIOD),
                                    Previous.size=seq(0,50,length=50),
                                    Reef.size = mean(dat$Reef.size)))
  mortframe$ctrsize <- mortframe$Previous.size-ctrval
  if (inherits(model,"mer")) {
    ## extract formula from model
    model.fixef.formula <- eval(lme4:::nobars(model@call$formula[-2]))
    X <- model.matrix(model.fixef.formula,data=mortframe)
    lmortvals <- X %*% fixef(model)
    ## predicted values
    mortvals <- plogis(lmortvals)
    var.fit <- diag(X %*% tcrossprod(vcov(model),X))
    se.fit <- sqrt(var.fit)
    var.tot <- var.fit + VarCorr(model)$REEF +VarCorr(model)$PR
    se.tot <- sqrt(var.tot)
    ## calc confidence intervals
    se.vals <- if (incl.var) se.tot else se.fit
    ## calc confidence intervals
    mort_lo <- plogis(lmortvals-1.96*se.vals)
    mort_hi <- plogis(lmortvals+1.96*se.vals)
  } else if (class(model)=="rjags") {
    X <- model.matrix(formula,data=mortframe)
    if (length(drop.cols)>0) X <- X[,-drop.cols]
    predmat <- plogis(apply(model$BUGSoutput$sims.list$beta,1,"%*%",t(X)))
    mortvals <- rowMeans(predmat)
    mort_lo <- apply(predmat,1,quantile,0.025)
    mort_hi <- apply(predmat,1,quantile,0.975)
  }
  mortframe <- data.frame(mortframe,mortvals,mort_lo,mort_hi)
  ##
  m.plot <-
    ggplot(dat,
           aes(x=Previous.size,
               y=Mortality,group=Treatment,
               colour=Treatment))+
                 geom_point(position=
                            position_jitter(height=jitter.y))+
                              facet_grid(~PERIOD)+
                                coord_cartesian(ylim=c(-exp.y,
                                                  1+exp.y))+
                                                    theme_bw()+
                                                      xlab("Previous size (cm)")+
                                                        ylab("Mortality probability")+
                                                          scale_colour_manual(values=colors)+
                                                            scale_fill_manual(values=colors)+
                                                      theme(title=element_text(title))
  ##
  m.plot2 <- m.plot+geom_line(data=mortframe,aes(y=mortvals),lwd=lwd)+
    geom_ribbon(data=mortframe,aes(y=mortvals,ymin=mort_lo,ymax=mort_hi,
                  fill=Treatment),
                alpha=0.3,colour=NA)
  if (ctrline) m.plot2 <- m.plot2 + geom_vline(xintercept=ctrval,col="darkgray",lty=2)
  m.plot2
}


coef.rjags <- function(object, type=c("mean","median"), ...) {
  type <- match.arg(type)
  cc <- object$BUGSoutput[[type]]
  cc$deviance <- NULL
  unlist(cc)
}

confint.rjags <- function(object, type=c("quantile","HPDinterval"), ...) {
  type <- match.arg(type)
  if (type=="HPDinterval") stop("not implemented")
  cc <- object$BUGSoutput$summary[,c("2.5%","97.5%")]
  cc <- cc[rownames(cc)!="deviance",]
  cc
}

## generate starting vals/inits
bfun <- function(drop.terms=character(0),dmat=model@X,model=m.acr.full,
                 sigma.pr=0.01,sigma.reef=0.01) {
  drop.cols <- colnames(dmat) %in% drop.terms
  X2 <- dmat
  if (length(drop.cols)>0) dmat <- dmat[,!drop.cols]
  b.dat <- with(model@frame,
               list(X2=X2,obs=Mortality,
                    pr=as.numeric(PR),reef=as.numeric(REEF),
                    n.pr=length(levels(PR)),
                    n.reef=length(levels(REEF)),
                    n.obs=length(Mortality),
                    n.beta=ncol(X2)))
  betavec <- fixef(model)
  if (length(drop.cols)>0) betavec <- betavec[!drop.cols]
  n.inits <- max(length(sigma.pr),length(sigma.reef))
  if (n.inits>1) {
    sigma.pr <- rep(sigma.pr,length=n.inits)
    sigma.reef <- rep(sigma.reef,length=n.inits)
  }
  b.inits <- mapply(function(p,r) {
    list(beta=betavec,sigma.pr=p,sigma.reef=r)
  },as.list(sigma.pr),as.list(sigma.reef),SIMPLIFY=FALSE)
  list(data=b.dat,inits=b.inits)
}

rename_coefs <- function(x,n,w) {
  x$BUGSoutput$mean <- as.list(unlist(x$BUGSoutput$mean)) ## could be dangerous?
  names(x$BUGSoutput$mean)[w] <- n
  colnames(x$BUGSoutput$sims.matrix)[w] <- n
  dimnames(x$BUGSoutput$sims.array)[[3]][w] <- n
  rownames(x$BUGSoutput$summary)[w] <- n
  x
}

histogram.mcmc <- function (x, data = NULL, outer,
                            aspect = "xy", default.scales = list(relation = "free"), 
                            start = 1, thin = 1, main = attr(x, "title"), xlab = "", 
                            ..., subset = coda:::thinned.indices(x, start = start, 
                                                        thin = thin)) 
{
    if (!is.R()) {
        stop("This function is not yet available in S-PLUS")
    }
    if (!missing(outer)) 
        warning("specification of outer ignored")
    data <- as.data.frame(x)
    form <- as.formula(paste("~", paste(lapply(names(data), as.name), 
        collapse = "+")))
    histogram(form, data = data[subset, , drop = FALSE], outer = TRUE, 
        aspect = aspect, default.scales = default.scales, main = main, 
        xlab = xlab, ...)
}

admbfun <- function(bdat,...) {
  tt <- system.time(d1 <- with(bdat,do_admb("demog1",data=data,
                                            params=inits[[1]],
                                            re=TRUE,
                                            re_vectors=c(u_rf=data$nrf,u_pr=data$npr),
                                            bounds=list("sigma_pr"=c(1e-4,20),"sigma_rf"=c(1e-4,20)),
                                            extra.args="-mno 1149",
                                            checkdata="write",
                                            checkparam="write",
                                            mcmcpars=c("beta","sigma_pr","sigma_rf"),
                                            mcmc=TRUE,...)))
  attr(d1,"runtime") <- tt
  d1
}
