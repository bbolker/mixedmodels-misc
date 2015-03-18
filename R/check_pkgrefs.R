library("tools")
library("stringr")

x <- readLines("../MixedModels.ctv")
a1 <- available.packages()

## package_dependencies defaults to Depends/Imports/LinkingTo
## could extend to "Suggests" as well, but the list is already pretty
##  darned long ...
pp <- unlist(package_dependencies(c("lme4","nlme","MCMCglmm"),
                           db=a1,
                           reverse=TRUE,recursive=FALSE))

##' extract one section from the text
##' @param tag XML section tag
##' @param x text
get_section <- function(tag,x) {
    beg <- grep(paste0("<",tag,">"),x)
    end <- grep(paste0("</",tag,">"),x)
    x[(beg+1):(end-1)]
}

##' remove beginning/ending whitespace, XML tags, and package delimiter ::
strip <- function(x) {
    gsub("(^ +| +$|<[^>]+>|::)","",x)
}

## 
pkglist <- strip(get_section("packagelist",x))
infolist <- strip(unlist(str_extract_all(get_section("info",x),
                              "<pkg>[^<]+</pkg>")))

cat("Packages in reverse depends list *not* in packagelist:\n")
cat(setdiff(pp,pkglist),"\n")
cat("Packages in packagelist *not* in reverse depends list:\n")
cat(setdiff(pkglist,pp),"\n")
cat("Packages in reverse depends list *not* in info:\n")
cat(setdiff(pp,infolist),"\n")
cat("Packages in info *not* in reverse depends list:\n")
cat(setdiff(infolist,pp),"\n")


                
## check
