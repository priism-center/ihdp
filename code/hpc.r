library(tidyverse)
library(rstan)
library(rstanarm)
library(survey)
source("code/matching.R")
source("code/balance.R")
source("code/estimation.R")
library(cobalt)
library(rgenoud)
library(Matching)
library(CBPS)
library(WeightIt)
library(foreach)
library(parallel)

load("data/cc2.Rdata")

covs_all <- setdiff(names(cc2), c("row.names", "row.names.1", "treat",
    "treat0", "ppvtr.36"))
covs <- c("bw", "preterm", "dayskidh", "sex", "first", "age", "black",
    "hispanic", "white", "b.marr", "lths", "hs", "ltcoll", "college",
    "work.dur", "prenatal", "momage")
cov_names <- c("birth weight", "weeks preterm", "days in hospital",
    "male", "first born", "age", "black", "hispanic", "white",
    "unmarried at birth", "less than high school", "high school graduate",
    "some college", "college graduate", "worked during pregnancy",
    "had no prenatal care", "age at birth")

ps_form <- formula(treat ~ bwg + hispanic + black + b.marr + lths + hs + ltcoll + work.dur + prenatal + sex + first + st5 + st9 + st12 + st25 + st36 + st42 + st48 + st53 + bw + preterm + momage + dayskidh)


set.seed(20)
        
cl <- makePSOCKcluster(10)
X <- setdiff(all.vars(ps_form), 'treat')

mgen_1_wr <- GenMatch(Tr=cc2$treat, X=cc2[,X], BalanceMatrix=cc2[,X], pop.size=1000, cluster=cl, balance=TRUE, replace=TRUE, ties=FALSE)

saveRDS(mgen_1_wr, model)
