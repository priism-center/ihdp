############################################
# Test script

############################################

library(rstanarm)
    cpus_per_task <- Sys.getenv('SLURM_CPUS_PER_TASK')
    if(cpus_per_task!=''){
        options(mc.cores=as.integer(cpus_per_task))
    } else {
        options(mc.cores=parallel::detectCores() - 1)
    }

############################################

set.seed(88)

df <- data.frame(x=rnorm(100, 8, 1), y=rnorm(100, 2*x, 3))

mod_test <- stan_glm(y ~ x, data=df, family=gaussian(link='identity'), algorithm='optimizing')

summary(mod_test)


############################################
# ps to te function

pscore_te <- function(spec, state=FALSE, te_spec=te_spec_nr){
    if (state){
        ps_spec <- update(spec, . ~ . + st5 + st9 + st12 + st25 + st36 + st42 + st48 + st53)
        te_spec <- update(te_spec, . ~ . + st5 + st9 + st12 + st25 + st36 + st42 + st48 + st53)
    }

    # Propensity score
    set.seed(8)
    ps_fit <- stan_glm(ps_spec, family=binomial(link='logit'), data=cc2, algorithm='optimizing')

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
