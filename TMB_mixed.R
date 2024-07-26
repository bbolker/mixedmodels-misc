library(TMB)
library(lme4)
library(glmmTMB)

## fit basic model 
m1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, REML = FALSE)
glmmtmb1 <- glmmTMB(Reaction ~ Days + (Days|Subject), sleepstudy, REML = FALSE)

## compile and load TMB DLL
compile("TMB_mixed.cpp") ##  "-O0 -g" for debugging
dyn.load(dynlib("TMB_mixed"))

## extract X and Z from lmer fit (lots of other ways to do this;
## model.matrix() is fine for X, we could use lme4::lFormula() for
## Z-construction instead of going all the way through lmer()
## if we want to build fancy RE model matrices ourselves we need
## Matrix::fac2sparse() and KhatriRao (see vignette("lmer", package = "lme4")
## for details)

X <- getME(m1, "X")
Z <- getME(m1, "Z")

## construct data and starting parameter values for TMB
tmbdat <- lme4:::namedList(X,
                           Z,
                           yobs = sleepstudy$Reaction,
                           n_re = 2L)

tmbpars <- list(beta = rep(0, ncol(X)),
                b = rep(0, ncol(Z)),
                theta = rep(0,3), ## 2 SD pars + 1 corr par ((n+1)*n/2)
                logsd = 0)

## build TMB object
obj <- MakeADFun(data = tmbdat,
                 parameters = tmbpars,
                 random = "b",
                 DLL = "TMB_mixed",
                 silent = TRUE ## FALSE for debugging etc.
                 )

## fit
tmbfit1 <- with(obj, nlminb(start = par, objective = fn, gradient = gr))

## checking glmmTMB vs lmer
## compare FE vcov
all.equal(as.matrix(vcov(m1)), vcov(glmmtmb1)$cond, tolerance = 1e-4)
## compare RE cov matrix
all.equal(c(VarCorr(m1)$Subject), c(VarCorr(glmmtmb1)$cond$Subject), tol = 1e-4)

all.equal(unname(glmmtmb1$fit$par),
          ## reorder parameters, and double logsd
          ##  (glmmTMB fits on the log-variance rather than the log-sd scale)
          unname(c(tmbfit1$par[1:2],
                   tmbfit1$par["logsd"]*2,
                   tmbfit1$par[3:5])),
          tolerance = 1e-6
          )

## try with BFGS
## ?? clearly suboptimal fit (5 log-likelihood units worse than tmbfit1 ...)
## would better starting values help??
## don't know what's going on here.
tmbfit2 <- with(obj, optim(par = par, fn = fn, gr = gr, method = "BFGS",
                           control = list(maxit = 200)))
