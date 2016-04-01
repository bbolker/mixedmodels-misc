---
Thoughts about design of statistical tools: what did R get right?
---

The [recent Moore foundation grant to Julia-lang](https://www.moore.org/newsroom/in-the-news/2015/11/10/bringing-julia-from-beta-to-1.0-to-support-data-intensive-scientific-computing) got me thinking about 

## update() with data by reference

It's quite common to fit many versions of the same model. If full or transformed versions of the data are stored inside the fitted model object (e.g. as model frames, design matrices, etc.), then this can quickly become unwieldy. It would be useful to a have a clean, universal design that allows groups of models to share a reference to the same data set. (This could be handled through a smart `update` method ...)

##  model structures: model matrices/frames/terms etc.

I've complained in the past about the confusion of model frames; `model.frame()` subsets to just the input variables that are used in the formula, omits records containing `NA`s in any of these columns. (So far so good; except in the case of fittine methods like random forests, or models with sets of regularized parameters like mixed effect models, we can't do anything with records with missing predictors anyway.) After that, it gets a little more mixed up. (1) Model terms involving data-dependent bases (splines, orthogonal polynomials, etc.) get computed from the data, and a `predvars` attribute gets attached to the `terms` attribute of the model frame to store the original information needed to recompute the basis. (2) Model terms with functions (e.g. `log(x)`) get evaluated. (3) Model terms involving factors do *not* get expanded; these remain as categorical variables until one calls `model.matrix`. I've asked Martin Maechler about the design decisions behind model matrices and terms, but I don't remember whether he gave me an answer.

All of this stuff is written in C for performance; in Julia it would presumably be written in Julia, which would be at least somewhat helpful for transparency and accessibility.

## contrasts

Perhaps contrasts are necessarily confusing, but R offers too many, non-transparent ways to change them. (1) global options; (2) as an attribute of a factor variable; (3) on the fly via `contrasts` in the model call (which gets passed to `model.matrix`, I guess ...) The naming convention for contrasts is confusing. (`car` offers some slight improvements.)

## make formulas first-class, manipulable objects

The "right" way to manipulate formulas is as language objects, rather than strings, but there aren't any helper functions to do so. `reformulate()` operates at a string level. `nlme`, `lme4`, and now `glmmTMB` have independently implemented, rather opaque recursive functions for transforming functions.

## flexible formulas

There are certain model structures that are very difficult to create via `model.matrix`. For example, suppose you set up a 2x2 factorial design ({`A`, `B`} x {`I`, `II`}) where you really want to remove the difference between I and II in the "A" treatment (suppose A and B are successive treatment periods, and the two treatments are randomized/normalized in the first period so their effect is 0 by definition). I don't know of any way to do this without creating a model matrix by hand and droppping columns.

There are also many, incompatible extensions of the formula language. `pscl` uses pipes to separate the conditional and zero-inflation parts of the model; `lme4` uses them to separate random-effect terms from random-effect grouping variables.

It might be nice if the variables output by `model.matrix` were themselves legal variable names.

## evaluation in different environments

- the R idiom of passing most parameters to a subsidiary function by changing the head of the function and evaluating it in the parent environment is sometimes fragile

## simulate methods

- better convention for storing "extra" parameters for simulation (e.g. NB or Gamma shape parameters, which are stored in very odd places)

## effects()/lsmeans equivalent

a convention for simulation/prediction of marginalizing and conditioning in different ways

## methods

- It would be nice if there were a uniform convention for parameter transformations (scaling, link/inverse-link transformations, orthogonalization)
- a convention for allowing the user to specify/fix some parameters ("masks")
- ditto, for `coef()`, a convention about subsets/classes of parameters for the same model (conditional vs zero-inflation parameters, fixed vs random vs conditional modes vs correlation/heteroscedasticity parameters, location vs scale parameters, etc)
- `simulate` and `predict` methods should allow new parameters to be specified (`newparams`) as well as new predictor-variable data (`newdata`); this enables methods like parametric bootstrapping

##  model types that can be separately defined and fitted

It's very useful to be able to set up an unevaluated model, optimizing it as a separate stage, as well as having a way to evaluate the objective (likelihood, posterior density, etc.) function of a model for user-specified parameters

## misc

[https://github.com/b-k/apophenia](Ben Klemens's apophenia modeling package) probably has some useful ideas
