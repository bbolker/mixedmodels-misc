---
name: MixedModels
topic: Mixed, multilevel, and hierarchical models in R
maintainer: Ben Bolker, ?
e-mail: bolker@mcmaster.ca, ?
version: 2022-07-29
source: https://github.com/bbolker/mixedmodels-misc/blob/master/taskview/MixedModels.md
---

**Authors**: Ben Bolker, Michael Agronah, ??

*Mixed models* are a broad class of statistical models used to analyze data where observations can be assigned to discrete groups, and where the parameters describing the differences are treated as *random variables*. They are also variously described as *multilevel*, *hierarchical*, or *repeated measures* models; *longitudinal* data are often analyzed in this framework as well.  Mixed models can be fitted in either frequentist or Bayesian frameworks.

**Scope**: only including models that incorporate *continuous* (usually although not always Gaussian) latent variables; this excludes packages that handle hidden Markov Models, finite (discrete) mixture models, latent Markov models, etc..

## Basic model fitting

LMMs are models with Normal residuals, responses that are linear combinations of the predictor variables, and Normal distributions of the random effects/latent variables.

### Linear mixed models

#### Frequentist

The most commonly used packages (functions) for frequentist LMMs are `r pkg("nlme", priority = "core")` (`lme`: REML/ML estimation, multiple nested random effects, residual correlation and heteroscedasticity) and `r pkg("lme4", priority = "core")` (`lmer`: REML/ML, nested and crossed REs, profile confidence intervals, parametric bootstrapping).

#### Bayesian 

Most packages use Markov chain Monte Carlo estimation: `r pkg("MCMCglmm", priority = "core")`, `r pkg("rstanarm")` and `r pkg("brms")`; the latter two packages are built on [Stan](mc-stan.org). `r pkg("blme")`, built on `r pkg("lme4")`, uses MAP estimation.
    
### Generalized linear mixed models

Generalized linear mixed models (GLMMs) can be seen as hierarchical extensions of generalized linear models (GLMs), or a extensions of LMMs to different response distributions (typically in the exponential family). The random-effect distributions are typically assumed to be Gaussian on the scale of the
linear predictor. 

#### Frequentist 

`r pkg("MASS")` (`glmmPQL`: fits via penalized quasi-likelihood), `r pkg("lme4")` (`glmer`: Laplace approximation and adaptive Gauss-Hermite quadrature), `r pkg("glmmTMB")` (Laplace approximation), `r pkg(GLMMadaptive)`, `r pkg("hglm")` (hierarchical GLMs). 
  
#### Bayesian 
       
`r pkg("MCMCglmm")`,`r pkg("rstanarm")`, `r pkg("brms")` and `r pkg("glmm")`  fit GLMMs in a Bayesian MCMC framework. `r pkg("MCMCglmm")`  fits GLMMs using Markov chain Monte Carlo techniques. `r pkg("rstanarm")` fits GLMMs using Markov Chain Monte Carlo, variational approximations to the posterior distribution, or optimization. `r pkg("brms")` supports a wide range of distributions and link functions for fitting GLMMs. `r pkg("glmm")` fits GLMMs using Monte Carlo Likelihood Approximation. 

**DELETE?** Packages specialized on binary data: `r pkg("glmmEP")` (expectation propagation, probit models only); `r pkg("GLMMRR")` (binary randomized response data)

### Nonlinear mixed models ([G]NLMMs)

NLMMs incorporate arbitrary nonlinear responses that cannot be accommodated in the framework of GLMMs. Only a few packages can accommodate **generalized** nonlinear mixed models
(i.e., nonlinear mixed models with non-Gaussian responses).

#### Frequentist

`r pkg("nlme")` (`nlme`), `r pkg("lme4")` (`nlmer`), `r pkg("repeated")` [GNLMM], `r pkg("saemix")` (stochastic approximation EM algorithm).

#### Bayesian

- `r pkg("brms")` [GNLMM?]

### Generalized estimating equations (GEEs)

GEEs represent an alternative approach to fitting clustered, longitudinal, or otherwise correlated data; GEE fits produce estimates of the *marginal* effects (averaged across the group-level variation) rather than *conditional* effects (conditioned on group-level information) `r pkg(" wgeesel")`, `r pkg("geesmv")`, `r pkg("geepack", priority = "core")`, `r pkg("gee")`, `r pkg("multgee")`  and `r pkg("geeM")`

## Specialized models

- **Robust estimation** (downweighting the importance of extreme observations: `r pkg("robustlmm")`, `r pkg("robustBLME")` (Bayesian robust LME), `r pkg("CRTgeeDR")`
- **Penalized models** (regularization or variable selection by ridge/lasso/elastic net penalties): `r pkg("splmm")` fits LMMs for high-dimensional data by imposing penalty on both the fixed effects and random effects for variable selection.
- **Handling missing values**: the `r pkg("mice")` package can be used to generate multiple imputation sets for use with other packages. `r pkg("mlmmm")` (EM imputation),  `r pkg("CRTgeeDR")` (GEEs )
- **Censored data** (responses : `r pkg("brms")` (general), `r pkg("lmec")` (censored Gaussian), `r pkg("ARpLMEC")` (censored Gaussian, autoregressive errors) and `r pkg("tlmec")` (censored Gaussian and t)
- **Ordinal-valued responses**: `r pkg("ordinal")`, `r pkg("cplm")`
- **Zero-inflated models**: (frequentist) `r pkg("glmmTMB")`, `r pkg("cplm")`; (Bayesian): `r pkg("MCMCglmm")`, `r pkg("brms")`.
- **Quantile regression**: `r pkg("lqmm")`, `r pkg("qrLM")`,`r pkg("qrNLMM")`
- **Phylogenetic/pedigree-based models**: `r pkg("pedigreemm")`, `r pkg("coxme")`, `r pkg("pez")`, `r pkg("kinship2")`
- **Survival analysis** (random effects are often referred to *frailty terms* in survival-analysis contexts): `r pkg("coxme")`
- **Spatial models**: [INLA](http://www.r-inla.org/home), `r pkg("nlme")` (with `corStruct` functions), `r pkg("CARBayesST")`, `r pkg("sphet")`, `r pkg("spind")`, `r pkg("spaMM")`, `r pkg("glmmfields")`, `r pkg("glmmTMB")`, `r pkg("inlabru")` (spatial point processes via log-Gaussian Cox processes) (brms?) (See also `r view("Handling and Analyzing Spatio-Temporal Data")`)
- **Differential equations**: `r pkg("mixedsde")`, `r pkg("nlmeODE")` `r pkg("PSM")`: see also `r view("Differential Equations")`
- **Large data sets**: `r pkg("mgcv")` (`bam()`)
- **Multinomial responses**: FIXME
- **Longitudinal data** (FIXME, explain): `r pkg("lmeNB")`
- **Factor analytic, latent variable, and structural equation modelling**:  `r pkg("lavaan")`, `r pkg("nlmm")`,`r pkg("sem")`, `r pkg("piecewiseSEM")`, `r pkg("semtree")`, `r pkg("semPLS")` and  `r pkg("blavaan")` . (See also the `r view("Psychometrics")` task view)
- **Tree-based models**: `r pkg("glmertree")`, `r pkg("semtree")`

## Model diagnostics and summary statistics

### Model diagnostics

`r pkg("HLMdiag")`, `r pkg("rockchalk")`, `r pkg("influence.ME")`, `r pkg("aods3")` (overdispersion), `r pkg("DHARMa")`, `r pkg("performance")` 

### Summary statistics

`r pkg("iccbeta")` (intraclass correlation), `r pkg("r2glmm")` (R^2 and partial R^2),
`r pkg("HiLMM")` (heritability), `r pkg("cAIC4")` (conditional AIC) , `r pkg("blmeco")` (WAIC)

### Derivatives

`r pkg("lmeInfo")`, `r pkg("merDeriv")`, `r pkg("lmmpar")`

**robust variance-covariance estimates**: `r pkg("clubSandwich")`, `r pkg("merDeriv")`

## Datasets

`r pkg("mlmRev")`, `r pkg("lme4")`, `r pkg("nlme")`, `r pkg("SASmixed")`, `r pkg("StroupGLMM")`, `r pkg("blmeco")`, `r pkg("nlmeU")`, `r pkg("VetResearchLMM")` 


## Model presentation and prediction

Functions and frameworks for convenient and tabular and graphical output of mixed model results: `r pkg("effects")` `r pkg("emmeans")`, `r pkg("dotwhisker")`, `r pkg("huxtable")`, `r pkg("sjPlots")`, `r pkg("rockchalk")`


`r pkg("broom.mixed")`, `r pkg("insight")`
 
## Convenience wrappers

These functions don't necessarily add new functionality, but 
provide convenient frameworks for less experienced users to fit and interpret mixed models.

 `r pkg("ez")`, `r pkg("mixlm")` `r pkg("afex")`, `r pkg("RVAideMemoire")`,
`r pkg("ZeligMultilevel")` `r pkg("cubature")`.


### Model and variable selection

`r pkg("LMERConvenienceFunctions")`, `r pkg("MuMIn")`

## Inference

`r pkg("pbkrtest")`, `r pkg("afex")` `r pkg("varTestnlme")`, `r pkg("lmeVarComp")`, `r pkg("RLRLsim")`, `r pkg("car")` (`Anova()`), `r pkg("CLME")`, `r pkg("lmertest")`

## Bootstrapping

`r pkg("pbkrtest")`, `r pkg("lme4")` (`bootMer` function), `r pkg("lmeresampler")`, `r pkg("glmmboot")`

#### Additive models

`r pkg("gamm4")`, `r pkg("mgcv")`, `r pkg("brms")`, `r pkg("lmeSplines")`

(note `assist` package is currently archived on CRAN)


### Bioinformatic applications

(FIXME: refer to/check bioconductor?)

`r pkg("MCMC.qpcr")`,`r pkg("CpGassoc")`, `r pkg("QGglmm")`, `r pkg("Phxnlme")`, `r pkg("mlmm.gwas")`

#### Power analysis

`r pkg("longpower")`, `r pkg("clusterPower")`, `r pkg("powerlmm")`, `r pkg("pass.lme")`.
     
#### Other

The following are other packages applied in mixed models. ; `r pkg("lmeNBBayes")` `r pkg("MarginalMediation")`
`r pkg("skewlmm")`  fits scale mixture of skew-normal linear mixed models using  expectation-maximization (EM) 
`r pkg("mvglmmRank")` implements multivariate Generalized Linear Mixed Models for ranking sport teams
   
## Links

- [R-SIG-mixed-models mailing list](https://stat.ethz.ch/mailman/listinfo/r-sig-mixed-models) for discussion of mixed-model-related questions, course announcements, etc.. 
- [r+mixed-models tags on Stack Overflow](http://stackoverflow.com/questions/tagged/r+mixed-models)
- [Cross Validated](http://stats.stackexchange.com)
- [INLA](http://www.r-inla.org/home)

## To do

- more cleanup
- more description?
- mention AS-REML, R interface; section on off-CRAN/commercial alternatives? (INLA, AS-REML, Mplus) Graphical front-ends (JASP/Jamovi)?
- missing packages: sommer, ... 
- scripts to check for archiving, download rank, etc.? link to/cross-check with [mixed model comparison table](https://docs.google.com/spreadsheets/d/19itelYaVW0U0gtNtRfqh76ZGt1awlamNcJwT71u_5Uk/edit#gid=0) ?
