---
name: MixedModels
topic: Mixed, multilevel, and hierarchical models in R
maintainer: Ben Bolker, ?
e-mail: bolker@mcmaster.ca, ?
version: 2022-07-29
source: https://github.com/bbolker/mixedmodels-misc/blob/master/taskview/MixedModels.md
---

*Mixed models* are a broad class of statistical models used to analyze data where observations can be assigned to discrete groups, and where the parameters describing the differences are treated as *random variables*. They are also variously described as *multilevel*, *hierarchical*, or *repeated measures* models; *longitudinal* data are often analyzed in this framework as well.  Mixed models can be fitted in either frequentist or Bayesian frameworks.

## Basic model fitting

LMMs are models with Normal residuals, responses that are linear combinations of the predictor variables, and Normal distributions of the random effects/latent variables.

### Linear mixed models

#### Frequentist

The most commonly used packages (functions) for frequentist LMMs are `r pkg("nlme")` (`lme`: REML/ML estimation, multiple nested random effects, residual correlation and heteroscedasticity) and `r pkg("lme4")` (`lmer`: REML/ML, nested and crossed REs, profile confidence intervals, parametric bootstrapping).

#### Bayesian 

Most packages use Markov chain Monte Carlo estimation: `r pkg("MCMCglmm")` (DESCRIBE), `r pkg("rstanarm")` and `r pkg("brms")`; the latter two packages are built on [Stan](mc-stan.org). `r pkg("blme")`, built on `r pkg("lme4")`, uses MAP estimation.
    
### Generalized linear mixed models

Generalized linear mixed models (GLMMs) can be seen as hierarchical extensions of generalized linear models (GLMs), or a extensions of LMMs to different response distributions (typically in the exponential family). The random-effect distributions are typically assumed to be Gaussian on the scale of the
linear predictor. 

#### Frequentist 

`r pkg("MASS")` (`glmmPQL`: penalized quasi-likelihood), `r pkg("lme4")` (`glmer`: Laplace approximation and adaptive Gauss-Hermite quadrature), `r pkg("glmmTMB")` (Laplace approximation), `r pkg(GLMMadaptive)`, `r pkg("hglm")` (hierarchical GLMs). 
  
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

GEEs are extensions of GLMs that fit longitudinal categorical responses data that are correlated. GEEs can also be applied in fit longitudinal continuous measurements. `r pkg(" wgeesel ")`,  `r pkg("geesmv")`, `r pkg("geepack")`, `r pkg(" gee")`, `r pkg(" multgee ")`  and `r pkg(" geeM ")` fit GEEs.   `r pkg(" wgeesel ")` fits weighted generalized estimating equations (WGEE). 

`r pkg("geesmv")` fits GEEs with more recent modified variance estimators for improving the finite small-sample performance.

`r pkg("geepack")` estimates parameters of the mean, scale and correlation structures of GEEs and can handle clustered categorical responses.

`r pkg(" gee ")` fits GEEs to data. User needs to specify any offsets in the model formula. 

`r pkg(" multgee ")` solves GEEs for correlated nominal or ordinal multinomial responses. `r pkg(" geeM ")` estimates GEEs parameters in mean structures with possible correlation between the outcomes.

## Specialized models

- **Robust estimation** (downweighting the importance of extreme observations: `r pkg("robustlmm")`, `r pkg("robustBLME")` (Bayesian robust LME), `r pkg("CRTgeeDR")`
- **Penalized models** (regularization or variable selection by ridge/lasso/elastic net penalties): `r pkg("splmm")` fits LMMs for high-dimensional data by imposing penalty on both the fixed effects and random effects for variable selection.
- **Handling missing values**: the `r pkg("mice")` package can be used to generate multiple imputation sets for use with other packages. `r pkg("mlmmm")` (EM imputation),  `r pkg("CRTgeeDR")` (GEEs )
- **Censored data** (responses : `r pkg("brms")` (general), `r pkg("lmec")` (censored Gaussian), `r pkg("ARpLMEC")` (censored Gaussian, autoregressive errors) and `r pkg("tlmec")` (censored Gaussian and t)
- **Ordinal data**: `r pkg("ordinal")`, `r pkg("cplm")`
- **Zero-inflated models**: (frequentist) `r pkg("glmmTMB")`, `r pkg("cplm")`; (Bayesian): `r pkg("MCMCglmm")`, `r pkg("brms")`.
- **Quantile regression**: `r pkg("lqmm ")`, `r pkg("qrLM")`,`r pkg("qrNLMM")`
- **Phylogenetic/pedigree-based models**: `r pkg("pedigreemm")`, `r pkg("coxme")`, `r pkg("pez")`, `r pkg("kinship2")`
- **Survival analysis** (random effects are often referred to *frailty terms* in survival-analysis contexts): `r pkg("coxme")`
 
## Model diagnostics and summary statistics

(See also "Inference")

### Model diagnostics

`r pkg("HLMdiag ")`, `r pkg("rockchalk")`, `r pkg("influence.ME")`, `r pkg("aods3")` (overdispersion), `r pkg("DHARMa")`, `r pkg("performance")` 

### Summary statistics

`r pkg("iccbeta ")` (intraclass correlation), `r pkg("r2glmm")` (R^2 and partial R^2),
`r pkg("HiLMM")` (heritability), `r pkg("cAIC4")` (conditional AIC) 

## Fitting Mixed Models to large dataset
These packages are appropriting for fitting mixed models when working with very large datasets. `r pkg("splmm ")` fit linear mixed-effects models for high-dimensional data imposing penalty on both the random and fixed effects for varaible selection. `r pkg("mgcv")`::` bam </code> fits generalized additive model to very large data set.   

## Dataset, sampling and functions for mixed models

The following packages contain datasets, functions and sampling methods for mixed model: 

<ul> <li> `r pkg("Blmeco ")`, `r pkg("nlmeU")`, `r pkg("SASmixed ")`, `r pkg("asremlPlus ")` and `r pkg("StroupGLMM ")` contain datasets for mixed model analysis . `r pkg("StroupGLMM ")` also contains R Codes for Generalized Linear Mixed Models  <li> `r pkg("Blmeco ")` and `r pkg("nlmeU")` also contain functions for analysing mixed models.  <li> `r pkg("lmeresampler")` and `r pkg("glmmboot")`: Implements bootstrap methods for mixed models  <li> `r pkg("arm")` has functions for computing the balance statistics.  <li> `r pkg("emmeans")` estimates marginal means and contrasts for mixed models <li> `r pkg("lsmeans")` computes least square means for mixed effect models <li> `r pkg("VetResearchLMM ")` contains Codes and Datasets for Linear Mixed Models in Veterinary Research. 

</ul>

## Model presentation and prediction

The following packages contain functions for convenient and tabular and graphical output of mixed model results: `r pkg("effects")` `r pkg("emmeans")`, `r pkg("dotwhisker")`, `r pkg("huxtable")`, `r pkg("sjPlots")`, `r pkg("rockchalk")`,`r pkg("broom.mixed")`
 
## Convenience wrappers

These functions don't necessarily add new functionality, but 
provide convenient frameworks for less experienced users to fit and interpret mixed models.

 `r pkg("ez")`, `r pkg("mixlm")` `r pkg("afex")`, `r pkg("RVAideMemoire")`,
`r pkg("ZeligMultilevel")` `r pkg("cubature ")`.


### Model and variable selection

`r pkg("LMERConvenienceFunctions")`, `r pkg("MuMIn")`

## Inference


The following packages are used for statistical inferences 

<ul> <li> `r pkg("asremlPlus ")` and `r pkg("pbkrtest")` contains functions for model testing  

<li> `r pkg("varTestnlme ")` impliments likelihood ratio Test for testing zero variances of a subset of the random effects in a mixed model.`r pkg("lmeVarComp ")` tests zero variance components in linear mixed models and test additivity in nonparametric regression  

<li> `r pkg("afex")` and `r pkg("ez")` provides functions for analyzing and visualising factorial experiments.

<li>
`r pkg("spaMM")` performs inference for mixed-effect models, including generalized linear mixed models with spatial correlations.

<li>
`r pkg("RLRsim")` implements a rapid simulation-based exact likelihood ratio tests
for testing the presence of variance /nonparametric terms for
models fit with nlme::lme(),lme4::lmer(), lmeTest::lmer(), gamm4::gamm4(),
mgcv::gamm() and SemiPar::spm(). 

<li> `r pkg("CLME")` performs inferencing for linear mixed models that have order constraints on some or all fixed effects. 

<li> `r pkg("lmerTest")` provides p-values in type I, II or III anova and summary tables for lmer model fits   

<li> `r pkg("lmem.qtler ")` performs QTL mapping analysis for balanced and for multi-environment and multi-trait analysis using mixed models  

<li>
`r pkg("rlme")` contains functions for estimating rank-based fixed effects and predicts robust random effects in two- and three- level random effects nested models. 

<li> `r pkg("MM4LMM ")` performs inference of Linear Mixed Models Through MM Algorithm 

</ul>

## Extensions
 
### Spatial/temporal models


  Geostatistical models (i.e. explicitly incorporating a model for continuous decay of correlation
 with distance, either in residuals or on random effects/latent variables); models based on *a priori*
 weights (e.g. simultaneous, conditional autoregression models). The following packages are used for Geostatistical models: <a href="http://www.r-inla.org/home">INLA</a>; `r pkg("nlme")` with `corStruct</code>; `r pkg("CARBayesST")`; `r pkg("sphet ")`; `r pkg("spind ")`, `r pkg("spaMM")` and  `r pkg("glmmfields ")`.
  
[INLA](http://www.r-inla.org/home) implements approximate Bayesian inference for Latent Gaussian Models; `r pkg("CARBayesST")` implements spatio-temporal generalised linear mixed models for areal unit data;  `r pkg("sphet ")` implements Generalized Method of Moment estimation for spatial autoregressive models with and without Heteroscedasticity. `r pkg("spind ")` has functions for spatial methods based on generalized estimating equations (GEE) and wavelet-revised methods (WRM) and  functions for spatially corrected model accuracy measures.

`r pkg("spaMM")` performce inference for mixed-effect models, including generalized linear mixed models with spatial correlations.  `r pkg("glmmfields ")` implements Bayesian spatial and spatiotemporal models and allows for extreme spatial deviations through time. 
 
   


User can specify a spatial correlation structure within the `corStruct</code>  argument of `r pkg("nlme")`. These are the available spatial correlation one can specify in `corStruct</code>


- `corExp</code>: exponential spatial correlation  
- `corGaus</code>: Gaussian spatial correlation  
-  `corLin</code>: linear spatial correlation  
- `corRatio</code>: Rational quadratics spatial correlation   
- `corSpher</code>: spherical spatial correlation 
</ul>



(See also "Handling and Analyzing Spatio-Temporal Data Task View") 

## Differential equation models

 The following are packages used for differential equation models with mixed effects: `r pkg("mixedsde ")`, `r pkg("nlmeODE")` and `r pkg("PSM ")`. `r pkg("mixedsde ")` is used for Inference on stochastic differential equation models invovling one or two mixed effects. `r pkg("PSM ")` provides functions for fitting linear and non-linear mixed-effects models using stochastic differential equations (SDEs). `r pkg("nlmeODE")` combines the odesolve and nlme packages for mixed-effects modelling using differential equations.     `r pkg("MsdeParEst")`: performs parametric estimation in stochastic differential equations with random effects in the drift, or in the diffusion or both  
 
(See also "Differential Equations task view")  

## Phylogenetic/pedigree-based models

A phylogenetic models study relationships among biological species based on similarities and differences in their physical or genetic characteristics. Pedigree models on the other hand, study the inheritance of a trait or disease accross several generations of biological species. 
 
The following packages are used for Phylogenetic/pedigree-based modeling.
 . `r pkg("pedigreemm")` fits LMMs or GLMMs incorporating the effects of pedigrees, `r pkg("coxme")` fit cox proportional hazards models containing both fixed and random effects and `r pkg("pez")` fits GLMMs for binary and continuous phylogenetic data. 
 

 
(See also "Phylogenetics task view")  

#### Additive models

`r pkg("gamm4")`, `r pkg("mgcv")`, `r pkg("brms")`, `r pkg("lmeSplines")`

(note `assist` package is currently archived on CRAN)


### Bioinformatic applications

These packages are useful mixed model packages applied in Bioinformatic: `r pkg("MCMC.qpcr")`,`r pkg("CpGassoc")`, `r pkg("QGglmm ")` and  `r pkg("Phxnlme ")`. 

  `r pkg("MCMC.qpcr")` analyses Quantitative RT-PCR data using GLMMs based on lognormal-Poisson error distribution, fitted using MCMC. `r pkg("CpGassoc")` can handle mixed effects models with chip or batch entering the model as a random intercept.

`r pkg("QGglmm ")` estimates Quantitative Genetics Parameters from Generalised Linear Mixed Models. `r pkg("Phxnlme ")` runs Phoenix NLME and Perform Post-Processing
Calls 'Phoenix NLME' (non-linear mixed effects), a population modeling and simulation software, for pharmacokinetics and pharmacodynamics analyses. 

### Factor analytic, latent variable, and structural equation models

Factor analytic models study variabilities in observed, correlated variables using unobserved variables known as factors. Latent variable models relates a set of observable variables (called manifest/response variables) to a set of latent variables (unobserved variables that are assumed to influence the response variables). Structural equation models combine a mixture of mathematical models, computer algorithms, and statistical methods in fitting networks of constructs to data.  

 The following packages are applied in factor analytic, latent variable, and structural equation modelling:  `r pkg("lavaan")`, `r pkg("nlmm ")`,`r pkg("sem")`, `r pkg("piecewiseSEM ")`, `r pkg("semtree")`, `r pkg("semPLS")` and  `r pkg("blavaan")` . 

 
 
`r pkg("lavaan")` fits a variety of latent variable models, including confirmatory factor analysis, structural equation and latent growth curve models.  `r pkg("nlmm ")` fit linear mixed models based on convolutions of the generalized Laplace (GL) distribution. 
 
 `r pkg("sem")`  conatins functions for fitting general linear structural
equation models with observed and latent variables. `r pkg("semtree")` constructs decision trees and forests to Structural Equation Models (SEM).

`r pkg("piecewiseSEM ")` implements piecewise structural equation modeling from a single list of structural equations and handles non-linear, latent, and composite variables, standardized coefficients, query-based prediction and indirect effects.  `r pkg("semPLS")` fits structural equation models using partial least squares (PLS).

`r pkg("blavaan")` fits a variety of Bayesian latent variable models, including confirmatory factor analysis, structural equation models, and latent growth curve models. 


 (See also "Psychometrics task view") 

#### Power analysis

Power analysis is used to explore how large a sample size needs to be in order to get a reasonably precise parameter estimates. The following packages are used for calculating power and sample sizes: `r pkg("longpower")`, `r pkg("clusterPower")`, `r pkg("powerlmm ")` and `r pkg("pass.lme")`.

    `r pkg("longpower")` computes power and sample size for LMMs and GAMs of longitudinal data, `r pkg("clusterPower")` computes power for cluster randomized trials (CRTs) GLMMs, `r pkg("powerlmm ")` calculates power for effects in multilevel longitudinal studies with missing data and  `r pkg("pass.lme")` computes power and sample size for testing fixed effect coefficients in multilevel linear mixed effect models.    

<bmb>descriptions!  emmeans does not fit models. We can probably take out glmmADMB. biglmm might belong in a separate "big data"category?</bmb> <mag> done </mag>

     
<h4>Other</h4>
The following are other packages applied in mixed models. `r pkg("blme")`; `r pkg("LMest ")`; `r pkg("lmeNBBayes ")`;`r pkg("lmeInfo ")`; `r pkg("MarginalMediation ")`; `r pkg("lmmpar ")`; `r pkg("ordinal")`. `r pkg("lmeInfo ")` provides analytic derivatives and information matrices for fitted linear mixed effects models; The `r pkg("ordinal")` package impliments cumulative link mixed models for analying ordinal (ordered categorical data) 


   `r pkg("skewlmm ")`  fits scale mixture of skew-normal linear mixed models using  expectation-maximization (EM) 

  
`r pkg("mvglmmRank ")` implements multivariate Generalized Linear Mixed Models for ranking sport teams ; `r pkg("mlmm.gwas")` implements multi-locus mixed model (MLMM) for  Genome-Wide Association Study (an area of study that detects common traits in a population and associations between genetic variants.). `r pkg(" glmertree ")` costructs tree models using GLMM or LMM for recurssive tree splits. `r pkg("lmeNB ")` fits longitudinal  count variables with a negative binomial mixed-effect regression model and uses maximum likelihood methods to estimate the  fixed effect parameters.
   
## Links

- [R-SIG-mixed-models mailing list](https://stat.ethz.ch/mailman/listinfo/r-sig-mixed-models) for discussion of mixed-model-related questions, course announcements, etc.. 
- [r+mixed-models tags on Stack Overflow](http://stackoverflow.com/questions/tagged/r+mixed-models)
- [Cross Validated](http://stats.stackexchange.com)
