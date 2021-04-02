MixedModels.html: MixedModels.ctv
	Rscript -e "ctv::ctv2html('MixedModels.ctv')"
	sed -i -e 's#\.\{2\}/packages#https://cran.r-project.org/web/packages#' MixedModels.html
## s=substitute
## #=delimiter
## between first and second delimiter = pattern to match
## between second and third delimiter = replacement string

## see https://regexr.com/
## or in R; x <- readLines("MixedModels.html"); gsub("\\.{2}packages","https://cran.r-project.org/web/packages",x);
##          writeLines(x,"MixedModels.html")

