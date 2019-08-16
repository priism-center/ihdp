############################################
# Exploring all possible propensity score models

############################################

library(rstan)
library(parallel)
    options(mc.cores=parallel::detectCores())
library(rstanarm)
library(survey)
library(cobalt)
source('library/matching.R')
source('library/balance.R')
source('library/estimation.R')

############################################
# Data

cc2 <- readRDS('data/ihdp.rds')

covs.all <- setdiff(names(cc2), c('row.names', 'row.names.1', 'treat', 'treat0', 'ppvtr.36'))
covs <- c('bw', 'preterm', 'dayskidh', 'sex', 'first', 'age', 'black', 'hispanic', 'white', 'b.marr', 'lths', 'hs', 'ltcoll', 'college', 'work.dur', 'prenatal', 'momage')
cov_names <- c('birth weight', 'weeks preterm', 'days in hospital', 'male', 'first born', 'age', 'black', 'hispanic', 'white', 'unmarried at birth', 'less than high school', 'high school graduate', 'some college', 'college graduate', 'worked during pregnancy', 'had no prenatal care', 'age at birth')

covs.nr <- c('bwg', 'hispanic', 'black', 'b.marr', 'lths', 'hs', 'ltcoll', 'work.dur', 'prenatal', 'sex', 'first', 'bw', 'preterm', 'momage', 'dayskidh')
covs.nr.st <- c(covs.nr, 'st5', 'st9', 'st12', 'st25', 'st36', 'st42', 'st53')

############################################
