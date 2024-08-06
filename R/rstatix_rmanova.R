## playing with classic ANOVA designs, comparing rstatix, aov, lme, lmerTest, ...

## what is rstatix anova_test() doing??

library(rstatix)
data("ToothGrowth")
df <- ToothGrowth
df$id <- rep(1:10, 6) # Add individuals id (fake, for RMANOVA purposes)
a1 <- anova_test(df,
           within = c("supp","dose"),
           wid = "id",
           dv = "len")
print(a1)
## sets up a multivariate linear model, does MANOVA ... calls car::Anova.mlm

## undebug(anova_test)
## undebug(.anova_test)
## debug(car_anova)
## debug(factorial_design)

## not sure how to specify that 
summary(aov(len~supp*dose + Error(id), data = df))
## summary(aov(len~supp*dose + Error(id/(supp*dose)), data = df))

## https://people.math.ethz.ch/~meier/teaching/anova/block-designs.html
## says: (fixed effect of id??)
summary(aov(len~id + supp*dose, data = df))


library(lmerTest)
m1 <- lmer(len ~ supp*dose + (supp*dose|id), data = df)
## anova() built into lmerTest: Satterthwaite by default
anova(m1)
anova(m1, ddf="Kenward-Roger")

## uses Kenward-Roger if 'F' is requested
car::Anova(m1, test.statistic = "F")


## nlme can't handle the overly complex random effects design
## could fart around with increasing number of iterations etc ...
library(nlme)
try(lme(len ~ supp*dose, random = ~supp*dose|id, data = df))
m2 <- lme(len ~ supp*dose, random = ~supp|id, data = df)
df$sd <- with(df, interaction(supp, dose))
m3 <- lme(len ~ supp*dose, random = ~1|id/sd, data = df)
anova(m3)

## ddf
-anova(m2)

source("https://raw.githubusercontent.com/bbolker/mixedmodels-misc/master/R/mm_utils.R")
calcDenDF(~supp*dose, random = ~supp*dose|id, data = df)
