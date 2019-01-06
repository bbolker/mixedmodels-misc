## new phyloglmm_setup

phylo.to.Z <- function(r,stand=FALSE){
  ntip <- length(r$tip.label)
  Zid <- Matrix(0.0,ncol=length(r$edge.length),nrow=ntip)
  nodes <- (ntip+1):max(r$edge)
  root <- nodes[!(nodes %in% r$edge[,2])]
  for (i in 1:ntip){
    cn <- i  ## current node
    while (cn != root){
      ce <- which(r$edge[,2]==cn)   ## find current edge
      Zid[i,ce] <- 1   ## set Zid to 1
      cn <- r$edge[ce,1]            ## find previous node
    }
  }
  V <- vcv(r)
  # V <- V/max(V)
  sig <- exp(as.numeric(determinant(V)["modulus"])/ntip)
  # sig <- det(V)^(1/ntip)
  Z <- t(sqrt(r$edge.length) * t(Zid))
  if(stand){Z <- t(sqrt(r$edge.length/sig) * t(Zid))}
  rownames(Z) <- r$tip.label
  colnames(Z) <- 1:length(r$edge.length)
  return(Z)                                  
}

phylo_lmm <- function(formula,data,phylo,phylonm=NULL,phyloZ=NULL,control,REML,doFit){
  lmod <- lFormula(formula=formula,data = data,control=control, REML=REML,phylonm=phylonm, phyloZ=phyloZ)
  devfun <- do.call(mkLmerDevfun, lmod)
  opt <- optimizeLmer(devfun,control=control$optCtrl)
  mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr)
}

lFormula <- function (formula, data = NULL, REML = TRUE, subset, weights, 
          na.action, offset, contrasts = NULL, control = lmerControl(),phylonm,phyloZ,  
          ...) 
{
  control <- control$checkControl
  mf <- mc <- match.call()
  ignoreArgs <- c("start", "verbose", "devFunOnly", "control")
  l... <- list(...)
  l... <- l...[!names(l...) %in% ignoreArgs]
  do.call(lme4:::checkArgs, c(list("lmer"), l...))
  if (!is.null(list(...)[["family"]])) {
    mc[[1]] <- quote(lme4::glFormula)
    if (missing(control)) 
      mc[["control"]] <- glmerControl()
    return(eval(mc, parent.frame()))
  }
  cstr <- "check.formula.LHS"
  lme4:::checkCtrlLevels(cstr, control[[cstr]])
  denv <- lme4:::checkFormulaData(formula, data, checkLHS = control$check.formula.LHS == 
                             "stop")
  formula <- as.formula(formula, env = denv)
  lme4:::RHSForm(formula) <- expandDoubleVerts(lme4:::RHSForm(formula))
  mc$formula <- formula
  m <- match(c("data", "subset", "weights", "na.action", "offset"), 
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  fr.form <- subbars(formula)
  environment(fr.form) <- environment(formula)
  for (i in c("weights", "offset")) {
    if (!eval(bquote(missing(x = .(i))))) 
      assign(i, get(i, parent.frame()), environment(fr.form))
  }
  mf$formula <- fr.form
  fr <- eval(mf, parent.frame())
  fr <- factorize(fr.form, fr, char.only = TRUE)
  attr(fr, "formula") <- formula
  attr(fr, "offset") <- mf$offset
  n <- nrow(fr)
  reTrms <- mkReTrms(findbars(lme4:::RHSForm(formula)), fr, phylonm, phyloZ)
  wmsgNlev <- lme4:::checkNlevels(reTrms$flist, n = n, control)
  wmsgZdims <- lme4:::checkZdims(reTrms$Ztlist, n = n, control, allow.n = FALSE)
  if (anyNA(reTrms$Zt)) {
    stop("NA in Z (random-effects model matrix): ", "please use ", 
         shQuote("na.action='na.omit'"), " or ", shQuote("na.action='na.exclude'"))
  }
  wmsgZrank <- lme4:::checkZrank(reTrms$Zt, n = n, control, nonSmall = 1e+06)
  fixedform <- formula
  lme4:::RHSForm(fixedform) <- nobars(lme4:::RHSForm(fixedform))
  mf$formula <- fixedform
  fixedfr <- eval(mf, parent.frame())
  attr(attr(fr, "terms"), "predvars.fixed") <- attr(attr(fixedfr, 
                                                         "terms"), "predvars")
  ranform <- formula
  lme4:::RHSForm(ranform) <- subbars(lme4:::RHSForm(lme4:::reOnly(formula)))
  mf$formula <- ranform
  ranfr <- eval(mf, parent.frame())
  attr(attr(fr, "terms"), "predvars.random") <- attr(terms(ranfr), 
                                                     "predvars")
  X <- model.matrix(fixedform, fr, contrasts)
  if (is.null(rankX.chk <- control[["check.rankX"]])) 
    rankX.chk <- eval(formals(lmerControl)[["check.rankX"]])[[1]]
  X <- lme4:::chkRank.drop.cols(X, kind = rankX.chk, tol = 1e-07)
  if (is.null(scaleX.chk <- control[["check.scaleX"]])) 
    scaleX.chk <- eval(formals(lmerControl)[["check.scaleX"]])[[1]]
  X <- lme4:::checkScaleX(X, kind = scaleX.chk)
  list(fr = fr, X = X, reTrms = reTrms, REML = REML, formula = formula, 
       wmsgs = c(Nlev = wmsgNlev, Zdims = wmsgZdims, Zrank = wmsgZrank))
}


mkReTrms <- function(bars, fr, phylonm,phyloZ,drop.unused.levels = TRUE){
  if (!length(bars)) 
    stop("No random effects terms specified in formula", 
         call. = FALSE)
  stopifnot(is.list(bars), vapply(bars, is.language, NA), inherits(fr, 
                                                                   "data.frame"))
  names(bars) <- lme4:::barnames(bars)
  term.names <- vapply(bars, lme4:::safeDeparse, "")
  blist <- lapply(bars, mkBlist, fr, phylonm, phyloZ, drop.unused.levels)
  nl <- vapply(blist, `[[`, 0L, "nl")
  if (any(diff(nl) > 0)) {
    ord <- rev(order(nl))
    blist <- blist[ord]
    nl <- nl[ord]
    term.names <- term.names[ord]
  }
  Ztlist <- lapply(blist, `[[`, "sm")
  Zt <- do.call(rBind, Ztlist)
  names(Ztlist) <- term.names
  q <- nrow(Zt)
  cnms <- lapply(blist, `[[`, "cnms")
  nc <- lengths(cnms)
  nth <- as.integer((nc * (nc + 1))/2)
  nb <- nc * nl
#   if (sum(nb) != q) {
#     stop(sprintf("total number of RE (%d) not equal to nrow(Zt) (%d)", 
#                  sum(nb), q))
#   }
  boff <- cumsum(c(0L, nb))
  thoff <- cumsum(c(0L, nth))
  Lambdat <- t(do.call(sparseMatrix, do.call(rBind, lapply(seq_along(blist), 
                                                           function(i) {
                                                             mm <- matrix(seq_len(nb[i]), ncol = nc[i], byrow = TRUE)
                                                             dd <- diag(nc[i])
                                                             ltri <- lower.tri(dd, diag = TRUE)
                                                             ii <- row(dd)[ltri]
                                                             jj <- col(dd)[ltri]
                                                             data.frame(i = as.vector(mm[, ii]) + boff[i], j = as.vector(mm[, 
                                                                                                                            jj]) + boff[i], x = as.double(rep.int(seq_along(ii), 
                                                                                                                                                                  rep.int(nl[i], length(ii))) + thoff[i]))
                                                           }))))
  thet <- numeric(sum(nth))
  ll <- list(Zt = drop0(Zt), theta = thet, Lind = as.integer(Lambdat@x), 
             Gp = unname(c(0L, cumsum(nb))))
  ll$lower <- -Inf * (thet + 1)
  ll$lower[unique(diag(Lambdat))] <- 0
  ll$theta[] <- is.finite(ll$lower)
  Lambdat@x[] <- ll$theta[ll$Lind]
  ll$Lambdat <- Lambdat
  fl <- lapply(blist, `[[`, "ff")
  fnms <- names(fl)
  if (length(fnms) > length(ufn <- unique(fnms))) {
    fl <- fl[match(ufn, fnms)]
    asgn <- match(fnms, ufn)
  }
  else asgn <- seq_along(fl)
  names(fl) <- ufn
  attr(fl, "assign") <- asgn
  ll$flist <- fl
  ll$cnms <- cnms
  ll$Ztlist <- Ztlist
  ll
}

mkBlist <- function (x, frloc, phylonm,phyloZ, drop.unused.levels = TRUE) 
{
  frloc <- factorize(x, frloc)
  if (is.null(ff <- tryCatch(eval(substitute(lme4:::makeFac(fac), 
                                             list(fac = x[[3]])), frloc), error = function(e) NULL))) 
    stop("couldn't evaluate grouping factor ", deparse(x[[3]]), 
         " within model frame:", " try adding grouping factor to data ", 
         "frame explicitly if possible", call. = FALSE)
  if (all(is.na(ff))) 
    stop("Invalid grouping factor specification, ", deparse(x[[3]]), 
         call. = FALSE)
  if (drop.unused.levels) 
    ff <- factor(ff, exclude = NA)
  if(phylonm[1] %in% names(frloc)){
    phyloZ <- phyloZ[levels(frloc[,phylonm[1]]),]
  }
  nl <- length(levels(ff))
  mm <- model.matrix(eval(substitute(~foo, list(foo = x[[2]]))), 
                     frloc)
  sm <- fac2sparse(ff, to = "d", drop.unused.levels = drop.unused.levels)
  if(grepl(phylonm[1],x[3])){
    # nbranch <- ncol(phyloZ)
    nrep <- nrow(sm)/nrow(phyloZ)
#     bmat <- matrix(0,nrow=nbranch,ncol=length(levels(ff)))
#     colnames(bmat) <- levels(ff)
#     for(i in rownames(phyloZ)){
#       bmat[,which(grepl(i,colnames(bmat)))] <- phyloZ[i,]
#     }
    lkr <- 1
    rkr <- 1
    if(nrep > 1){
      if(strsplit(as.character(x[[3]]),":")[[2]] == phylonm[1]){
        lkr <- nrep
      }
      if(strsplit(as.character(x[[3]]),":")[[2]] == phylonm[2]){
        rkr <- nrep
      }
    # sm <- t(sm)
    }
    # bmat <- kronecker(kronecker(diag(lkr),phyloZ),diag(rkr))
    bmat <- kronecker(diag(lkr),phyloZ)
 #   if(nrow(sm) == nspp*nsite){
  #    sm <- t(sm[ff,])
  #  }
    sm <- t(t(sm) %*% bmat)
    nl <- nrow(sm)
  }
  sm <- KhatriRao(sm, t(mm))
  # dimnames(sm) <- list(rep(1:nrow(sm), each = ncol(mm)), rownames(mm))
  list(ff = ff, sm = sm, nl = nl, cnms = colnames(mm))
}


phylo_glmm <- function(formula,data,phylo,phylonm=NULL,phyloZ=NULL,control,family){
  glmod <- glFormula(formula=formula,data = data,control=control,family,phylonm=phylonm, phyloZ=phyloZ)
  # glmod$reTrms <- modify_phylo_retrms(glmod$reTrms,phylo,phylonm,phyloZ)
  devfun <- do.call(mkGlmerDevfun, glmod)
  opt <- optimizeGlmer(devfun)
  devfun <- updateGlmerDevfun(devfun,glmod$reTrms)
  opt <- optimizeGlmer(devfun,stage=2)
  mkMerMod(environment(devfun), opt, glmod$reTrms, fr = glmod$fr)
}

glFormula <- function (formula, data = NULL, family = gaussian, subset, weights, 
                        na.action, offset, contrasts = NULL, start, mustart, etastart, 
                        control = glmerControl(),phylonm,phyloZ,...) 
{
  control <- control$checkControl
  mf <- mc <- match.call()
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame(2))
  if (is.function(family)) 
    family <- family()
  if (isTRUE(all.equal(family, gaussian()))) {
    mc[[1]] <- quote(lme4::lFormula)
    mc["family"] <- NULL
    return(eval(mc, parent.frame()))
  }
  if (family$family %in% c("quasibinomial", "quasipoisson", 
                           "quasi")) 
    stop("\"quasi\" families cannot be used in glmer")
  ignoreArgs <- c("start", "verbose", "devFunOnly", "optimizer", 
                  "control", "nAGQ")
  l... <- list(...)
  l... <- l...[!names(l...) %in% ignoreArgs]
  do.call(lme4:::checkArgs, c(list("glmer"), l...))
  cstr <- "check.formula.LHS"
  lme4:::checkCtrlLevels(cstr, control[[cstr]])
  denv <- lme4:::checkFormulaData(formula, data, checkLHS = control$check.formula.LHS == 
                             "stop")
  mc$formula <- formula <- as.formula(formula, env = denv)
  m <- match(c("data", "subset", "weights", "na.action", "offset", 
               "mustart", "etastart"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  fr.form <- subbars(formula)
  environment(fr.form) <- environment(formula)
  for (i in c("weights", "offset")) {
    if (!eval(bquote(missing(x = .(i))))) 
      assign(i, get(i, parent.frame()), environment(fr.form))
  }
  mf$formula <- fr.form
  fr <- eval(mf, parent.frame())
  fr <- factorize(fr.form, fr, char.only = TRUE)
  attr(fr, "formula") <- formula
  attr(fr, "offset") <- mf$offset
  if (!missing(start) && is.list(start)) {
    attr(fr, "start") <- start$fixef
  }
  n <- nrow(fr)
  reTrms <- mkReTrms(findbars(lme4:::RHSForm(formula)), fr,phylonm, phyloZ)
  wmsgNlev <- lme4:::checkNlevels(reTrms$flist, n = n, control, allow.n = TRUE)
  wmsgZdims <- lme4:::checkZdims(reTrms$Ztlist, n = n, control, allow.n = TRUE)
  wmsgZrank <- lme4:::checkZrank(reTrms$Zt, n = n, control, nonSmall = 1e+06, 
                          allow.n = TRUE)
  fixedform <- formula
  lme4:::RHSForm(fixedform) <- nobars(lme4:::RHSForm(fixedform))
  mf$formula <- fixedform
  fixedfr <- eval(mf, parent.frame())
  attr(attr(fr, "terms"), "predvars.fixed") <- attr(attr(fixedfr, 
                                                         "terms"), "predvars")
  ranform <- formula
  lme4:::RHSForm(ranform) <- subbars(lme4:::RHSForm(lme4:::reOnly(formula)))
  mf$formula <- ranform
  ranfr <- eval(mf, parent.frame())
  attr(attr(fr, "terms"), "predvars.random") <- attr(terms(ranfr), 
                                                     "predvars")
  X <- model.matrix(fixedform, fr, contrasts)
  if (is.null(rankX.chk <- control[["check.rankX"]])) 
    rankX.chk <- eval(formals(lmerControl)[["check.rankX"]])[[1]]
  X <- lme4:::chkRank.drop.cols(X, kind = rankX.chk, tol = 1e-07)
  if (is.null(scaleX.chk <- control[["check.scaleX"]])) 
    scaleX.chk <- eval(formals(lmerControl)[["check.scaleX"]])[[1]]
  X <- lme4:::checkScaleX(X, kind = scaleX.chk)
  list(fr = fr, X = X, reTrms = reTrms, family = family, formula = formula, 
       wmsgs = c(Nlev = wmsgNlev, Zdims = wmsgZdims, Zrank = wmsgZrank))
}
