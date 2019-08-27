############################################
# Personal R library

############################################

library(parallel)
# Setup cores
if(Sys.getenv("SLURM_CPUS_PER_TASK") != "") {
    options(mc.cores=as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")))
} else {
  options(mc.cores=parallel::detectCores()-1)
}

library(cobalt)


############################################


#' Find the ps_spec which gave greatest balance
#' 
#' @param i numeric, size of covariate sets to search
#' @param mwr bool, TRUE use MWR sample
findMax <- function(i, mwr=TRUE){
    bals <- readRDS(paste0('outputs/ps_bals_i',i,'.rds'))
    specs <- readRDS(paste0('outputs/ps_specs_i',i,'.rds'))

    if (mwr){
        bals_m <- lapply(bals, function(x) x[2])
    } else {
        bals_m <- lapply(bals, function(x) x[1])
    }

    idx <- which.max(bals_m)
    return(specs[idx])
}


#' Specific IHDP balance table
#' 
#' @param data data.frame of matched sample
bal_ihdp <- function(data){
    # Covariates to check balance on
    covs <- c('bw', 'bwg', 'preterm', 'dayskidh', 'sex', 'first', 'black', 'hispanic', 'b.marr', 'lths', 'hs', 'ltcoll', 'work.dur', 'prenatal', 'momage', 'income')
    bal.tab(data=data, treat=data$treat, covs=data[,covs], estimand='ATT', m.threshold=0.1, v.threshold=1.1, continuous='std', binary='raw')
}


#' Propensity score matching treatment effect estimation
#' 
#' @param spec formula; propensity score specification
#' @param state lgl; use state indicators in both ps & te specs
#' @param te_spec formula; treatment effect specficiation
#' @return data.frame of MwoR & MwR TEs
pscore_te <- function(spec, state=FALSE, te_spec=te_spec_nr){
    if (state){
        spec <- update(spec, . ~ . + st5 + st9 + st12 + st25 + st36 + st42 + st48 + st53)
        te_spec <- update(te_spec, . ~ . + st5 + st9 + st12 + st25 + st36 + st42 + st48 + st53)
    }

    # Propensity score estimation
    set.seed(8)
    ps_fit <- stan_glm(spec, family=binomial(link='logit'), data=cc2, algorithm='optimizing')

    pscores <- apply(posterior_linpred(ps_fit, type='link'), 2, mean)
    
    # Matching
    matches <- matching(z=cc2$treat, score=pscores, replace=FALSE)
    matches_wr <- matching(z=cc2$treat, score=pscores, replace=TRUE)
    
    # Restructuring sample
    matched <- cc2[matches$match.ind,]
    matched_wr <- cc2[matches_wr$match.ind,]

    # Treatment effect estimation
    set.seed(8)
    reg_ps <- stan_glm(te_spec, data=matched, algorithm='optimizing')
    
    reg_ps_wr.design <- svydesign(ids=~1, weights=matches_wr$cnts, data=cc2)
    reg_ps_wr <- svyglm(te_spec, design=reg_ps_wr.design, data=cc2)

    # Estimates
    te <- summary(reg_ps)['treat', 1:2]
    te_wr <- summary(reg_ps_wr)$coef['treat', 1:2]

    as.data.frame(c(te, te_wr))
}
