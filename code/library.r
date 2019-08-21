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
