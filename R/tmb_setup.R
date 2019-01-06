### TMB Helper functions

modify_TMBstruc <- function(tmbstruc,phylo,phylonm,
                            phyloZ) {
  n.edge <- nrow(phylo$edge)
  ## stuff in tmbstruc (but not within data.tmb) is (maybe) necessary
  ##  for cosmetics, proper reporting
  tmbstruc$condList$reTrms <- modify_phylo_retrms(tmbstruc$condList$reTrms, #tmbstruc$condList$reTrms,
                                                  phylo,phylonm,phyloZ)
  tmbstruc$condReStruc$`1 + X | sp`$blockReps <- n.edge
  tmbstruc$condList$Z <- t(tmbstruc$condList$reTrms$Zt)
  ## data *inside* data.tmb is actually the most critical to allow correct fit
  tmbstruc$data.tmb$terms$`1 + X | sp`$blockReps <- n.edge
  tmbstruc$data.tmb$Z <- t(tmbstruc$condList$reTrms$Zt)
  tmbstruc$parameters$b <- rep(0,ncol(tmbstruc$data.tmb$Z))
  return(tmbstruc)
}

fit_TMBstruc <- function(TMBStruc,verbose=FALSE) {
  obj <- with(TMBStruc, TMB:::MakeADFun(data.tmb, parameters, map = mapArg, 
    random = randomArg, profile = NULL, silent = !verbose, DLL = "glmmTMB"))
    optTime <- system.time(fit <- with(obj, nlminb(start = par, objective = fn, 
                                                   gradient = gr)))
    fit$parfull <- obj$env$last.par.best
    sdr <- TMB:::sdreport(obj)
    report <- obj$report()
    return(list(fit=fit,sdr=sdr,report=report))
  }