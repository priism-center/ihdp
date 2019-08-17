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
# Functions


#' Estimate propensity score, run balance diagnostic, return aggregate balance metric
#' 
#' @param spec model specification
#' @param data data
#' @param covs covariates for balance assessment
#' @param mt mean threshold
#' @param vt variance threshold
#' @return numeric vector of number of total covariates meeting mean & variance thresholds in matched w/o replacement and matched w/ replacement
psBal <- function(spec, data=cc2, covs=covs_bal, mt=0.1, vt=1.1){
    set.seed(8)

    ps_fit <- stan_glm(spec, data=data, algorithm='optimizing')
    pscores <- apply(posterior_linpred(ps_fit, type='link'), 2, mean)
    matches <- matching(z=cc2$treat, score=pscores, replace=FALSE)
    matches_wr <- matching(z=cc2$treat, score=pscores, replace=TRUE)
    matched <- data[matches$match.ind,]
    matched_wr <- data[matches_wr$match.ind,]

    bal <- bal.tab(data=matched, treat=matched$treat, covs=matched[,covs], estimand='ATT', m.threshold=mt, v.threshold=vt, continuous='std', binary='raw')
    bal_wr <- bal.tab(data=matched_wr, treat=matched_wr$treat, covs=matched_wr[,covs], estimand='ATT', m.threshold=mt, v.threshold=vt)

    bal_metric <- bal$Balanced.Means$count[1] + bal$Balanced.Variances$count[1]
    bal_metric_wr <- bal_wr$Balanced.Means$count[1] + bal_wr$Balanced.Variances$count[1]

    return(c(bal_metric, bal_metric_wr))
}


############################################
# Data

cc2 <- readRDS('data/ihdp.rds')

cc2$bwT = (cc2$bw-1500)^2
cc2$dayskidT = log(cc2$dayskidh+1)
cc2$pretermT = (cc2$preterm+8)^2
cc2$momageT = (cc2$momage^2)

# Covariates to balance
covs_bal <- c('bw', 'bwg', 'preterm', 'dayskidh', 'sex', 'first', 'black', 'hispanic', 'b.marr', 'lths', 'hs', 'ltcoll', 'work.dur', 'prenatal', 'momage', 'income')
covs_bal_st <- c(covs_bal, 'st5', 'st9', 'st12', 'st25', 'st36', 'st42', 'st48', 'st53')

# Covariates for propensity score estimation
covs_ps <- c('bw', 'preterm', 'dayskidh', 'sex', 'first', 'bwg', 'black', 'hispanic', 'b.marr', 'lths', 'hs', 'ltcoll', 'work.dur', 'prenatal', 'momage', 'income', 'bwT', 'dayskidT', 'pretermT', 'momageT', 'income')
# Create all two way two-way interactions
n <- length(covs_ps)
idx <- unlist(lapply(1:n, function(i) combn(1:n, 2, simplify=FALSE)),recursive=FALSE)
covs_int <- sapply(idx, function(i) paste(covs_ps[i], collapse=':'))


############################################
# Create all propensity score specifications
covs_1 <- c(covs_ps, covs_int)

n <- length(covs_1)

idx <- unlist(
    lapply(1:n, function(i)
        combn(1:n, i, simplify=FALSE)
    ), recursive=FALSE)

if (file.exists('outputs/ps_specs.rds')){
    ps_specs <- readRDS('outputs/ps_specs.rds')
} else {
    ps_specs <- sapply(idx, function(i)
        paste0('treat~', paste(covs_1[i], collapse='+'))
        )
}

# Estimation
if (file.exists('outputs/ps_bals.rds')){
    ps_bals <- readRDS('outputs/ps_bals.rds')
} else {
    ps_bals <- mclapply(ps_specs, function(spec) psBal(spec), mc.cores=detectCores())
}

# Save everything
saveRDS(ps_specs, 'outputs/ps_specs.rds')
saveRDS(ps_bals, 'outputs/ps_bals.rds')
