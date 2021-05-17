## https://www.dropbox.com/s/0w3csm2vfcj47yq/data_email.csv?dl=0

library(MASS)
library(lme4)
library(bbmle)
library(sjPlot)
library(ggplot2); theme_set(theme_bw())
library(ggeffects)

## power-logistic link function
logexp <- function(expos = 1)  {
   linkfun <- function(mu) qlogis(mu^(1/expos))
   ## FIXME: is there some trick we can play here to allow
   ##   evaluation in the context of the 'data' argument?
   linkinv <- function(eta)  plogis(eta)^expos
   mu.eta <- function(eta) expos * plogis(eta)^(expos-1) *
     .Call(stats:::C_logit_mu_eta, eta, PACKAGE = "stats")
   valideta <- function(eta) TRUE
   link <- paste("logexp(", deparse(substitute(expos)), ")",
                 sep="")
   structure(list(linkfun = linkfun, linkinv = linkinv,
                  mu.eta = mu.eta, valideta = valideta,
                  name = link),
             class = "link-glm")
}

FL1 <- read.csv("logexp_data.csv")

## original goal (fails)
try(g250.3 <-glmer(fate~(1|brood)+log(I(1+G250))+(log(I(1+G250^2)))+(log(I(1+G250^3))),
               family=binomial(link=logexp(FL1$Interval)), data=FL1, start=NULL))

## check balance/sample size
with(FL1,table(fate))


## cloglog models (linear, polynomial)
gm1 <- glmer(fate~(1|brood) + G250 + offset(log(Interval)),
             family=binomial(link="cloglog"),
             data=FL1)
gm2 <- update(gm1, . ~ . -G250 + poly(G250,2))
gm3 <- update(gm1, . ~ . -G250 + poly(G250,3)) ## warning

## try with power-logistic
gm2B <- update(gm1, family=binomial(link=logexp(FL1$Interval)))
a1 <- allFit(gm2B) ## looks OK-ish
summary(a1)$fixef  ## coefficient comparisons

## quadratic model fails though (even with starting vals derived from lin model)
try(gm2C <- update(gm1,
               . ~ . -G250 + poly(G250,2),
               family=binomial(link=logexp(FL1$Interval)),
               start=list(beta=c(fixef(gm2B,0)))))

## compare what we've got so far
bbmle::AICctab(gm1, gm2, gm3, gm2B, mnames=c("linear","quad","cubic", "lin_logexp"))

## plot/compare predictions
g1 <- ggpredict(gm1, terms="G250 [0:8.5 by=0.2]")
## predict by hand for logexp, it's hard to update the exposure term in the link function
X <- cbind(int=1,G250=seq(0,8.5,by=0.2))
linkinv <- logexp(expos=mean(FL1$Interval))$linkinv
g2B <- data.frame(x=X[,"G250"], predicted=  linkinv(X %*% fixef(gm2B)))
g2B.sd <- sqrt(diag(X %*% vcov(gm2B) %*% t(X)))
g2B <- within(g2B, {
              conf.low <- predicted-1.96*g2B.sd
              conf.high <- predicted+1.96*g2B.sd
              })

## now plot results together
## (better to get all predictions, rbind/bind_rows together, plot from scratch)
## slightly cobbled together because of how sjPlot does things
sjPlot::plot_model(gm2, type="pred", terms="G250 [0:8.5 by=0.2]") +
  geom_line(data=g1, colour="red") +
  geom_ribbon(data=g1, colour=NA, fill="red", alpha=0.2,
              aes(ymin=conf.low, ymax=conf.high)) +
  geom_line(data=g2B, colour="blue") +
  geom_ribbon(data=g2B, colour=NA, fill="blue", alpha=0.2,
              aes(ymin=conf.low, ymax=conf.high)) +
  stat_sum(data=FL1, alpha=0.5, aes(x=G250, y=fate))

