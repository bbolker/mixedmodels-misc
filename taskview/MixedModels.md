---
name: MixedModels
topic: Mixed, multilevel, and hierarchical models in R
maintainer: Ben Bolker, ?
e-mail: bolker@mcmaster.ca, ?
version: 2022-07-29
---

*Mixed models* are a broad class of statistical models used to analyze data where observations can be assigned to discrete groups, and where the parameters describing the differences are treated as *random variables*. They are also variously described as *multilevel*, *hierarchical*, or *repeated measures* models; *longitudinal* data are often analyzed in this framework as well.  Mixed models can be fitted in either frequentist or Bayesian frameworks.

The [R-SIG-mixed-models mailing list](https://stat.ethz.ch/mailman/listinfo/r-sig-mixed-models) is an active forum for discussion of mixed-model-related questions, course announcements, etc.. [Stack Overflow](http://stackoverflow.com/questions/tagged/r+mixed-models) and [Cross Validated](http://stats.stackexchange.com) also host relevant discussions.

## Model fitting

LMMs are models with Normal residuals, responses that are linear combinations of the predictor variables, and Normal distributions of the random effects/latent variables.

### Linear mixed models

#### Frequentist

The most common packages (functions) for frequentist LMMs are `r pkg("nlme")` (`lme`: REML/ML estimation, multiple nested random effects, residual correlation and heteroscedasticity) and `r pkg("lme4")` (`lmer`: REML/ML, nested and crossed REs, profile confidence intervals, parametric bootstrapping).

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

# ## Generalized estimating equations
GEEs are extensions of GLMs that fit longitudinal categorical responses data that are correlated. GEEs can also be applied in fit longitudinal continuous measurements. `r pkg(" wgeesel ")`,  `r pkg("geesmv")`, `r pkg("geepack")`, `r pkg(" gee")`, `r pkg(" multgee ")`  and `r pkg(" geeM ")` fit GEEs.   `r pkg(" wgeesel ")` fits weighted generalized estimating equations (WGEE). 

`r pkg("geesmv")` fits GEEs with more recent modified variance estimators for improving the finite small-sample performance.

`r pkg("geepack")` estimates parameters of the mean, scale and correlation structures of GEEs and can handle clustered categorical responses.

`r pkg(" gee ")` fits GEEs to data. User needs to specify any offsets in the model formula. 


`r pkg(" multgee ")` solves GEEs for correlated nominal or ordinal multinomial responses. `r pkg(" geeM ")` estimates GEEs parameters in mean structures with possible correlation between the outcomes.

### Special topics

- **Robust estimation** (downweighting the importance of extreme observations: `r pkg("robustlmm")`, `r pkg("robustBLME")` (Bayesian robust LME), `r pkg("CRTgeeDR")`
- **Penalized models** (regularization or variable selection by ridge/lasso/elastic net penalties): `r pkg("splmm")` fits LMMs for high-dimensional data by imposing penalty on both the fixed effects and random effects for variable selection.
- **Handling missing values**: the `r pkg("mice")` package can be used to generate multiple imputation sets for use with other packages. `r pkg("mlmmm")` (EM imputation),  `r pkg("CRTgeeDR")` (GEEs )
- **Censored data** (responses : `r pkg("brms")` (general), `r pkg("lmec")` (censored Gaussian), `r pkg("ARpLMEC")` (censored Gaussian, autoregressive errors) and `r pkg("tlmec")` (censored Gaussian and t)
 
## Model diagnostics and summary statistics
The following packages provide tools and functions for models diagnotics
 <li>`r pkg("HLMdiag ")`, `r pkg("asremlPlus ")` and `r pkg("rockchalk ")` contains model diagnostics functions  <li> 
`r pkg("influence.ME ")` provides tools for identifying influential data points in mixed effect models 
<li> `r pkg("aods3 ")` is used for analysing overdispersed counts or proportions  <li> `r pkg("r2glmm")`: computes R squared for linear and generalized linear mixed effect models 
<li> `r pkg("iccbeta ")` used to compute intraclass correlation in grouped data 
<li> `r pkg("DHARMa ")` computes scaled residuals for fitted (generalized) linear mixed models for easiness of interpretation <li> `r pkg("qrLMM ")`: implements quantile Regression for Linear Mixed-Effects Models  
<li> `r pkg("HiLMM ")`: Estimates Heritability with confidence intervals in linear mixed models 
</ul>



(See also "Inference")
## Fitting Mixed Models to large dataset
These packages are appropriting for fitting mixed models when working with very large datasets. `r pkg("splmm ")` fit linear mixed-effects models for high-dimensional data imposing panalty on both the random and fixed effects for varaible selection. `r pkg("mgcv")`::` bam </code> fits generalized additive model to very large data set.   


 

## Dataset, sampling and functions for mixed models

The following packages contain datasets, functions and sampling methods for mixed model: 

<ul> <li> `r pkg("Blmeco ")`, `r pkg("nlmeU")`, `r pkg("SASmixed ")`, `r pkg("asremlPlus ")` and `r pkg("StroupGLMM ")` contain datasets for mixed model analysis . `r pkg("StroupGLMM ")` also contains R Codes for Generalized Linear Mixed Models  <li> `r pkg("Blmeco ")` and `r pkg("nlmeU")` also contain functions for analysing mixed models.  <li> `r pkg("lmeresampler")` and `r pkg("glmmboot")`: Implements bootstrap methods for mixed models  <li> `r pkg("arm")` has functions for computing the balance statistics.  <li> `r pkg("emmeans")` estimates marginal means and contrasts for mixed models <li> `r pkg("lsmeans")` computes least square means for mixed effect models <li> `r pkg("VetResearchLMM ")` contains Codes and Datasets for Linear Mixed Models in Veterinary Research. 

</ul>

## Model presentation and prediction
The following packages contain fuctions tabular and graphical representing of mixed model results as well as other forms. 

<ul> <li> `r pkg("effects")` creates graphical and tabular effect display for linear models  <li> `r pkg("dotwhisker")` quick and easy tool for ploting box and wiskers for model results  <li> `r pkg("huxtable")` creates tables for latex and HTML and other formats. `r pkg("sjPlots")`, `r pkg("rockchalk")` and `r pkg("asremlPlus ")` creates figures and tables for data visualization  <li> `r pkg("broom.mixed")` converts object from some of the mixed model packages in R into data frames  <li> `r pkg("car")` computes type-II or type-III analysis-of-variance tables  
</ul> 
 
 ## Convenience wrappers
These functions don't necessarily add new functionality, but 
provide convenient frameworks for (typically) less experienced
users to fit and interpret mixed models.
 `r pkg("ez")`, `r pkg("mixlm")` `r pkg("afex")`, `r pkg("RVAideMemoire")`,
`r pkg("ZeligMultilevel")` `r pkg("cubature ")`.  `r pkg("pez")` provide wrappers for common community phylogenetic indices.


### Model and variable selection

i> `r pkg("nmle::")``  anova.lme </code> and  `r pkg("nmle::")`` anova.gls </code> functions compare the likelihoods of multple fitted models.  
  
  <li>
`r pkg("ASReml-R ")` provides information criteria for selecting terms to include in models fitted by `r pkg("asreml ")`  
<li> `r pkg("cAIC4 ")`: Estimates conditional Akaike information for generalized mixed-effect models 
<li> stepAIC in the Mass package 
`r pkg("LMERConvenienceFunctions ")` is used for model selection  <li> `r pkg("MuMIn ")`Provides information criteria for model selection and carry out model averaging based on information criteria  <li> `r pkg("glmmLasso ")` is used for variable selection for generalized linear mixed models using L1-Penalization. FlexParamCurve for model section.  
<li> ` nls_multstart </code> function in the `r pkg("nls.multstart ")` package
finds the best fit of non-linear model based on AIC score  
 </ul>
 
##  Quantile Regression Mixed Models
 Quantile regression models estimate the conditional median (quantiles) of the response accross the values of the covariates. The following packages are used for fitting Quantile Regression Models involving fixed and random  effects: `r pkg("lqmm ")` and  `r pkg("qrNLMM ")`. 
 
 
`r pkg("lqmm ")` fits quantile regression models for hierarchical data. 
 `r pkg("qrNLMM ")` preforms Quantile regression (QR) for Nonlinear Mixed-Effects Models.    

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
  
 
  
 <a href="http://www.r-inla.org/home">INLA</a> implements approximate Bayesian inference for Latent Gaussian Models; `r pkg("CARBayesST")` implements spatio-temporal generalised linear mixed models for areal unit data;  `r pkg("sphet ")` implements Generalized Method of Moment estimation for spatial autoregressive models with and without Heteroscedasticity. `r pkg("spind ")` has functions for spatial methods based on generalized estimating equations (GEE) and wavelet-revised methods (WRM) and  functions for spatially corrected model accuracy measures.

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
`r pkg("pedigreemm")`, `r pkg("coxme")` and `r pkg("pez")`  . `r pkg("pedigreemm")` fits LMMs or GLMMs incorporating the effects of pedigrees, `r pkg("coxme")` fit cox proportional hazards models containing both fixed and random effects and `r pkg("pez")` fits GLMMs for binary and continuous phylogenetic data. 
 

 
(See also "Phylogenetics task view")  

<h4>Generalized additive models, splines</h4>
Generalized additive models (GAMs) are extensions of generalized linear model where the linear response variable depends linearly on unknown smooth functions to be estimated in the model. The following packages are used for fitting GAMs to data: `r pkg("gamm4")` (depends on `r pkg("lme4 (>= 1.0)")`, `r pkg("mgcv (>= 1.7-23)) ")`, `r pkg("gam")`, `r pkg("mgcv")``::gamm</code> (depends on `r pkg("nlme ((>= 3.1-64)")`) and `r pkg("gamlss ")`.   

 

 `r pkg("gamm4")`  is based on from package `r pkg("mgcv ")`::` gamm </code> and fits GAMs to data using fitting functions from `r pkg("lme4 ")`.  Selection of the Smoothness function is done by  Restricted Maximum Likelihood (REML) in the Gaussian additive case and (Laplace approximate) Maximum Likelihood otherwise. 

`r pkg("mgcv")`::` gam </code> fits GAMs with integrated smoothness estimation and implements  Generalized Cross Validation  for choosing the degree of freedom of the smoothing functions.
`r pkg("mgcv")`::` bam </code> fits generalized additive model to very large data set. 



`r pkg("lmeSplines")` fits smoothing spline terms in Gaussian linear and nonlinear mixed-effects models; the ` slm </code> and ` snm </code> functions in the `r pkg("assist ")` package fit semiparametric linear mixed-effects models and semiparametric nonlinear mixed-effects models respectively.  



<h4>Bioinformatic applications</h4>
These packages are useful mixed model packages applied in Bioinformatic: `r pkg("MCMC.qpcr")`,`r pkg("CpGassoc")`, `r pkg("QGglmm ")` and  `r pkg("Phxnlme ")`. 

  `r pkg("MCMC.qpcr")` analyses Quantitative RT-PCR data using GLMMs based on lognormal-Poisson error distribution, fitted using MCMC. `r pkg("CpGassoc")` can handle mixed effects models with chip or batch entering the model as a random intercept.

`r pkg("QGglmm ")` estimates Quantitative Genetics Parameters from Generalised Linear Mixed Models. `r pkg("Phxnlme ")` runs Phoenix NLME and Perform Post-Processing
Calls 'Phoenix NLME' (non-linear mixed effects), a population modeling and simulation software, for pharmacokinetics and pharmacodynamics analyses. 


<url> https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html </url>

## Zero-inflated models

The following packages are used for zero-inflated models (models with  probability distributions that allows for frequent zero-valued observations): <a href="https://glmmadmb.r-forge.r-project.org/">glmmADMB</a>, `r pkg("MCMCglmm")` and `r pkg("cplm")`.

The <a href="https://glmmadmb.r-forge.r-project.org/">glmmADMB</a>  package can handle simple (intercept-only) zero-inflated mixed models.`r pkg("MCMCglmm")` handles zero-inflated, hurdle, and zero-altered models (with arbitrary models for zeros) by stacking observations and constructing a multi-type model (see sections 5.3-5.5 of the `CourseNotes</code> vignette). The `r pkg("cplm")` package provides both frequentist and Bayesian (MCMC) tools for fitting zero-inflated compound Poisson (Tweedie) mixed models.


## Ecological and environmental applications

These packages are applied in ecological and environmental modeling

 `r pkg("HydroME")`::`SSomuto</code> ; `r pkg("dsm")`;   
 `r pkg("carcass")`::` search.efficiency</code>; `r pkg("blmeco")` and `r pkg("secr")`
 
 
 `r pkg("HydroME")`::`SSomuto</code> fits water retention characteristics for a grouped dataset as well as with mixed-effects modelling. 
 
 `r pkg("carcass")`::` search.efficiency</code>  estimates detection probabilty of the number of fatalities from carcass searches using a binomial model with vegetation density as fixed effect and person as random factor. 
 
`r pkg("blmeco")`::` WAIC ` computes Widely Applicable Information criterion (WAIC) for measuring predictive fit for mixed models.
 
 
`r pkg("dsm")` fits a density surface model (DSM) to detection adjusted counts from a spatially-referenced distance sampling analysis using generalized mixed models or generalized additive models.

 `r pkg("secr")` contains functions to estimate the density and size of a spatially distributed animal population. 


   (See also "Environmental task view") 
 
 <h4>Factor analytic, latent variable, and structural equation models</h4>
Factor analytic models study variabilities in observed, correlated variables using unobserved variables known as factors. Latent variable models relates a set of observable variables (called manifest/response variables) to a set of latent variables (unobserved variables that are assumed to influence the response variables). Structural equation models combine a mixture of mathematical models, computer algorithms, and statistical methods in fitting networks of constructs to data.  

 The following packages are applied in factor analytic, latent variable, and structural equation modelling:  `r pkg("lavaan")`, `r pkg("nlmm ")`,`r pkg("sem")`, `r pkg("piecewiseSEM ")`, `r pkg("semtree")`, `r pkg("semPLS")` and  `r pkg("blavaan")` . 

 
 
`r pkg("lavaan")` fits a variety of latent variable models, including confirmatory factor analysis, structural equation and latent growth curve models.  `r pkg("nlmm ")` fit linear mixed models based on convolutions of the generalized Laplace (GL) distribution. 
 
 `r pkg("sem")`  conatins functions for fitting general linear structural
equation models with observed and latent variables. `r pkg("semtree")` constructs decision trees and forests to Structural Equation Models (SEM).

`r pkg("piecewiseSEM ")` implements piecewise structural equation modeling from a single list of structural equations and handles non-linear, latent, and composite variables, standardized coefficients, query-based prediction and indirect effects.  `r pkg("semPLS")` fits structural equation models using partial least squares (PLS).

`r pkg("blavaan")` fits a variety of Bayesian latent variable models, including confirmatory factor analysis, structural equation models, and latent growth curve models. 


 (See also "Psychometrics task view") 
<h4>Power analysis</h4>
Power analysis is used to explore how large a sample size needs to be in order to get a reasonably precise parameter estimates. The following packages are used for calculating power and sample sizes: `r pkg("longpower")`, `r pkg("clusterPower")`, `r pkg("powerlmm ")` and `r pkg("pass.lme")`.

    `r pkg("longpower")` computes power and sample size for LMMs and GAMs of longitudinal data, `r pkg("clusterPower")` computes power for cluster randomized trials (CRTs) GLMMs, `r pkg("powerlmm ")` calculates power for effects in multilevel longitudinal studies with missing data and  `r pkg("pass.lme")` computes power and sample size for testing fixed effect coefficients in multilevel linear mixed effect models.    

<bmb>descriptions!  emmeans does not fit models. We can probably take out glmmADMB. biglmm might belong in a separate "big data"category?</bmb> <mag> done </mag>

     
<h4>Other</h4>
The following are other packages applied in mixed models. `r pkg("blme")`; `r pkg("LMest ")`; `r pkg("lmeNBBayes ")`;`r pkg("lmeInfo ")`; `r pkg("MarginalMediation ")`; `r pkg("lmmpar ")`; `r pkg("ordinal")`. `r pkg("lmeInfo ")` provides analytic derivatives and information matrices for fitted linear mixed effects models; The `r pkg("ordinal")` package impliments cumulative link mixed models for analying ordinal (ordered categorical data) 


   `r pkg("skewlmm ")`  fits scale mixture of skew-normal linear mixed models using  expectation-maximization (EM) 

  
`r pkg("mvglmmRank ")` implements multivariate Generalized Linear Mixed Models for ranking sport teams ; `r pkg("mlmm.gwas")` implements multi-locus mixed model (MLMM) for  Genome-Wide Association Study (an area of study that detects common traits in a population and associations between genetic variants.). `r pkg(" glmertree ")` costructs tree models using GLMM or LMM for recurssive tree splits. `r pkg("lmeNB ")` fits longitudinal  count variables with a negative binomial mixed-effect regression model and uses maximum likelihood methods to estimate the  fixed effect parameters.
   
 
  

