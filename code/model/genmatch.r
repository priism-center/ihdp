#!/share/apps/r/3.6.0/intel/bin/R

############################################
# Genetic matching on IHDP data

############################################

args <- commandArgs(trailingOnly=TRUE)

if (length(args)==0){
    REP <- FALSE
    model_file <- 'models/mgen_1.rds'
} else {
    REP <- TRUE
    model_file <- 'models/mgen_1_wr.rds'
}

############################################

library(parallel)

############################################

load('data/cc2.Rdata')

ps_form <- formula(treat ~ bwg + hispanic + black + b.marr + lths + hs + ltcoll + work.dur + prenatal + sex + first + st5 + st9 + st12 + st25 + st36 + st42 + st48 + st53 + bw + preterm + momage + dayskidh)

############################################

set.seed(20)
cl <- makePSOCKcluster(detectCores(), type='PSOCK')

X <- setdiff(all.vars(ps_form), 'treat')

mgen_mod <- GenMatch(Tr=cc2$treat, X=cc2[,X], BalanceMatrix=cc2[,X], pop.size=1000, cluster=cl, replace=REP)

saveRDS(mgen_mod, model_file)
