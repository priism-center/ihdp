############################################
# Propensity score matched treatment effect estimation

############################################

library(rstan)
    options(mc.cores=parallel::detectCores())
library(arm)
library(rstanarm)
library(survey)
library(parallel)
source('library/matching.R')
source('library/balance.R')
source('library/estimation.R')

############################################

load('data/cc2.Rdata')


covs.all <- setdiff(names(cc2), c('row.names', 'row.names.1', 'treat', 'treat0', 'ppvtr.36'))
covs <- c('bw', 'preterm', 'dayskidh', 'sex', 'first', 'age', 'black', 'hispanic', 'white', 'b.marr', 'lths', 'hs', 'ltcoll', 'college', 'work.dur', 'prenatal', 'momage')
cov_names <- c('birth weight', 'weeks preterm', 'days in hospital', 'male', 'first born', 'age', 'black', 'hispanic', 'white', 'unmarried at birth', 'less than high school', 'high school graduate', 'some college', 'college graduate', 'worked during pregnancy', 'had no prenatal care', 'age at birth')

covs.nr <- c('bwg', 'hispanic', 'black', 'b.marr', 'lths', 'hs', 'ltcoll', 'work.dur', 'prenatal', 'sex', 'first', 'bw', 'preterm', 'momage', 'dayskidh')
covs.nr.st <- c(covs.nr, 'st5', 'st9', 'st12', 'st25', 'st36', 'st42', 'st53')

############################################

psEst2 <- function(te_spec, dataMWOR=matched, matchesWR=matches_wr, data=cc2){
    set.seed(8)
    te_reg <- stan_glm(te_spec, data=matched, algorithm='optimizing')
    te_reg_wr_design <- svydesign(ids=~1, weights=matches_wr$cnts, data=data)
    te_reg_wr <- svyglm(te_spec, design=te_reg_wr_design, data=data)

    mwor <- summary(te_reg)['treat', 1:2]
    mwr <- summary(te_reg_wr)$coef['treat', 1:2]

    rbind(mwor, mwr)
}

############################################

ps_spec <- formula(treat ~ bwg*as.factor(educ) + as.factor(ethnic)*b.marr + work.dur + prenatal + preterm + momage + sex + first + bw + dayskidT + pretermT + momageT + black*(bw + preterm + dayskidT) + b.marr*(bw + preterm + dayskidT))

set.seed(8)

    ps_fit <- stan_glm(ps_spec, family=binomial(link='logit'), data=cc2, algorithm='optimizing')
    pscores <- apply(posterior_linpred(ps_fit, type='link'), 2, mean)
    matches <- matching(z=cc2$treat, score=pscores, replace=FALSE)
    matches_wr <- matching(z=cc2$treat, score=pscores, replace=TRUE)
    matched <- cc2[matches$match.ind,]
    matched_wr <- cc2[matches_wr$match.ind,]

############################################
n <- length(covs.nr)

# All possible combinations
id <- unlist(
        lapply(1:n, function(i)
            combn(1:n, i, simplify=FALSE)
        )
      , recursive=FALSE)

# Formulas
te_specs <- sapply(id, function(i)
              paste("ppvtr.36 ~ treat + ",paste(covs.nr[i],collapse="+"))
            )

# estimation
numCores <- detectCores()
ests <- mclapply(te_specs, function(spec) psEst2(spec), mc.cores=numCores)
