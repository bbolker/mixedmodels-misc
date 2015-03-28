library("tools")
library("stringr")
library("sos")
library("XML")
library("ctv")

if (FALSE) {
    mm_ctv <- read.ctv("MixedModels.ctv")
    ctv2html(mm_ctv)
    browseURL("./MixedModels.html")
}

mm_text <- readLines("../MixedModels.ctv")
a1 <- available.packages()

## package_dependencies defaults to Depends/Imports/LinkingTo
## could extend to "Suggests" as well, but the list is already pretty
##  darned long ...
pp <- unlist(package_dependencies(c("lme4","nlme","MCMCglmm"),
                           db=a1,
                           reverse=TRUE,recursive=FALSE))


dd <- findFn("glmm")
class(dd) <- "data.frame"  ## strip/avoid browser window
sos_pkgs <- unique(dd$Package)
pp <- union(pp,sos_pkgs)



theurl <- "http://cran.r-project.org/web/packages/available_packages_by_name.html"
descr_table_full <- readHTMLTable(theurl,colClasses=rep("character",2))[[1]]
descr_table <- setNames(descr_table_full[descr_table_full[[1]] %in% pp,],
                        c("pkg","description"))
## can View() or write.csv()
    
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
pkglist <- strip(get_section("packagelist",mm_text))
infolist <- strip(unlist(str_extract_all(get_section("info",x),
                              "<pkg>[^<]+</pkg>")))

cat("Packages in scraped package list *not* in packagelist:\n")
cat(setdiff(pp,pkglist),"\n")
cat("Packages in packagelist *not* in scraped package list:\n")
cat(setdiff(pkglist,pp),"\n")
cat("Packages in scraped package list *not* in info:\n")
cat(setdiff(pp,infolist),"\n")
cat("Packages in info *not* in scraped package list:\n")
cat(setdiff(infolist,pp),"\n")
cat("Packages in info *not* in packagelist:\n")
cat(setdiff(infolist,pkglist),"\n")

                
## check

