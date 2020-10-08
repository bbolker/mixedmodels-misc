## utilities/misc code for finding interesting packages related to mixed models

library(tidyverse) ## general 
library(miniCRAN)
library(crandep)
library(igraph)
## library(packdep) ## archived ...
library(packageRank)

## 1. plot reverse dependencies of lme4
## dependency types to include
## A depends on B: B automatically gets loaded
## A imports B: A uses functions from B, B *must* be installed
## A suggests B: A optionally uses functions from B
## A enhances B: the authors of the package think these go together
## "reverse-depends": reverse dependencies of lme4 = all the packages that depend on lme4
rd <- c("Reverse depends", "Reverse imports", "Reverse linking to", "Reverse suggests")
dd <- get_dep_df("lme4",rd)
plot(igraph::graph_from_data_frame(dd))  ## too many!
## grep/regular expression tricks

## 2. find many (not all) GLMM packages
a1 <- available.packages()
## grep("lmm",rownames(a1),value=TRUE,ignore.case=TRUE)
## lm followed by ("m or e") followed by (a character not "t" or the end of the string)
regexps <- c("lm(m|e([^t]|$))")  ## was using "mixed" but ... ? should check,
find_pkgs <- function(x) grep(x,rownames(a1),value=TRUE,ignore.case=TRUE)
focal_pkgs <- character(0)
for (r in regexps) {
    focal_pkgs <- union(focal_pkgs, find_pkgs(r))
}
## false pos
fpos <- "palmerpenguins"
## false negatives: some known-interesting pkgs
## (check MixedModels.ctv for some more)
fneg <- c("SASmixed","broom.mixed",
          "pbkrtest","emmeans","mgcv","gamm4",
          "brms","rstanarm","pez","merDeriv","repeated","hglm",
          "geesmv","geepack","influence.ME","cAIC4","HLMdiag","lmmfit","iccbeta",
          "DHARMa","effects","rockchalk","arm","performance","car",
          "ez","afex","RVAideMemoire","geoRglm","GLMMarp","spaMM",
          "polytomous","ordinal","longpower")
focal_pkgs <- focal_pkgs %>% setdiff(fpos) %>% union(fneg)
length(focal_pkgs) ## 112

## now extract
pkg_rd <- (expand_grid(name=focal_pkgs, type=rd))

## clunky
pb <- txtProgressBar(max=nrow(pkg_rd),style=3)
i <- 0
ff <- function(name,type) {
    ## cat(".")
    cat(name,"\n")
    i <<- i+1
    ## setTxtProgressBar(pb,i)
    tibble(focal=name,type,dep_pkg=get_dep(name,type))
}


if (file.exists("all_deps.rds")) {
    all_deps <- readRDS("all_deps.rds")
} else {
    ## THIS BIT IS SLOW, WATCH OUT ... (~ 6 minutes)
    system.time(all_deps <- (pkg_rd
        %>% pmap(ff)
        %>% bind_rows()
    )
    )
    saveRDS(all_deps,"all_deps.rds")
}

## don't drop NA vals yet -- about 2/3 of packages (all but 28) have no RDs at all

unique_deps <- (all_deps
    %>% drop_na(dep_pkg)
    %>% select(focal,dep_pkg)
    %>% unique()
)

ud2 <- (unique_deps
    %>% filter(dep_pkg %in% focal_pkgs)
)

rdg1 <-igraph::graph_from_data_frame(ud2[,c(2,1)])
plot(rdg1)
     
## 3. collect importance measures

cc <- eigen_centrality(rdg1)$vector
central_tbl <- tibble(focal=names(cc),central=cc)

a1a <- (all_deps
    %>% mutate_at("type",str_remove,"Reverse ")
    %>% mutate_at("type",
                  ~ case_when(. %in% c("depends","imports") ~ "strong",
                              TRUE ~ "weak"))
)

a1 <- (a1a
    %>% drop_na()
    %>% count(focal, type)
)

a_tot <- a1  %>% group_by(focal) %>% summarise(n=sum(n),.groups="drop") %>% mutate(type="total")

a2 <- (bind_rows(a1,a_tot)
    %>% pivot_wider(names_from=type, values_from=n,values_fill=0)
    ## restore packages with no depends
    %>% full_join(tibble(focal=focal_pkgs),by="focal")
    %>% mutate_all(replace_na,0)
)

## SLOW the first time
pp <- packageRank(focal_pkgs)
pp2 <- (pp$package.data
    %>% as_tibble()
    %>% select(focal=packages,downloads,percentile)
)

a3 <- (full_join(pp2, a2, by="focal")
    %>% full_join(central_tbl, by="focal")
    %>% mutate_at("central",replace_na,0)
    %>% mutate(score=percentile/100+strong/max(strong)+central)
    %>% arrange(desc(score))
)

a3
View(a3)
## add descriptions??

write_csv(a3,"glmm_packages.csv")




