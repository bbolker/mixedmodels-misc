## There are newer versions of packages replacing 'packdep':
## crandep
## pkggraph
if (FALSE) {
    library(dlstats)
    library(ggplot2); theme_set(theme_bw())
    cc <- cran_stats(c("crandep","pkggraph"))
    ggplot(cc, aes(start,downloads,colour=package)) + geom_line()
}    

library(pkggraph)
library(igraph)

## general-purpose "apply across combinations of multiple lists"
xapply <- function(FUN,...,FLATTEN=TRUE,MoreArgs=NULL) {
  ## add progress bar??
  L <- list(...)
  inds <- do.call(expand.grid,lapply(L,seq_along)) ## Marek's suggestion
  retlist <- vector("list",nrow(inds))
  for (i in 1:nrow(inds)) {
    arglist <- mapply(function(x,j) x[[j]],L,as.list(inds[i,]),SIMPLIFY=FALSE)
    if (FLATTEN) {
      retlist[[i]] <- do.call(FUN,c(arglist,MoreArgs))
    }
  }
  retlist
}

pkggraph::init(local=FALSE,repos="CRAN")
g1 <- get_neighborhood("nlme",)

tmpf <- function(repos,pkg,order) {
    ## repos ignored: could be Bioconductor/CRAN
    g <- make_neighborhood_graph(get_neighborhood(pkg,level=order))
    n <- igraph::ego_size(g[[1]],match(pkg,V(g[[1]])$name),order=order,mode="in")
    return(n)
}
dd <- expand.grid(repos=c("CRAN"),pkg=c("nlme","lme4"),order=1:3)
dd$number <- unlist(xapply(tmpf,
                           list(dd),list("nlme","lme4"),as.list(1:3)))

library(ggplot2); theme_set(theme_bw())
ggplot(dd,aes(order,number,colour=pkg,lty=repos,shape=repos))+
    geom_point()+geom_line()+
        theme_bw()+scale_x_continuous(breaks=1:3)+
            scale_y_log10()

png("lme4dep1.png",1200,1200)
gg <- make_neighborhood_graph(get_neighborhood("lme4",level=1))
plot(gg)
dev.off()

