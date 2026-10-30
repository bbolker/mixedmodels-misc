#' implementation of "between-within" df-guessing algorithm
#' based on Pinheiro and Bates 2000 description, but tries
#' to do better with random-slopes models

#' @examples
#' calcDenDF(~age,"Subject",nlme::Orthodont)
#' calcDenDF(~age,data=nlme::Orthodont,random=~1|Subject)
#' calcDenDF(~age,data=nlme::Orthodont,random=~age|Subject)
#' ## off by 1
#' calcDenDF(~poly(age,2),data=nlme::Orthodont,
#'            random=~poly(age,2)|Subject)
#' library(nlme)
#' lmeDF <- function(formula=distance~age,random=~1|Subject) {
#'     mod <- lme(formula,random,data=Orthodont)
#'     aa <- anova(mod)
#'     return(setNames(aa[,"denDF"],rownames(aa)))
#' }
#' lmeDF()
#' lmeDF(random=~age|Subject) ## wrong!
#' @param fixed fixed-effect formula
#' @param grps (character) vector of grouping vectors
#' @param random random-effect formula (can be specified instead of grps)
#' @param n0 ?? number of parameters fitted at the population level
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
