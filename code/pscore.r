############################################
# Functions to abstract and simplify the pscore modeling process
############################################

library(rstan)
    options(mc.cores=parallel.detectCores())
libary(rstanarm)
library(arm)
source('code/matching.R')
source('code/balance.R')
source('code/estimationn.R')

############################################

load('data/cc2.Rdata')


# Create pscore, do matching, create balance plot, estimate treatment effect on restructured sample
pscoreAll <- function(covs, form=NULL, matching=FALSE, data=cc2){

    if (is.null(form)) {
        ps_form <- as.formula(data[, c('treat', covs)])
    } else {
        ps_form <- form
    }
    ps_fit <- stan_glm(ps_form, family=binomial(link='logit'), data=data, algorithm='optimizing')

    pscores <- apply(posterior_linpred(ps_fit, type='link'), 2, mean)
    matches <- matching(z=data$treat, score=pscores, replace=matching)

    bal <- balance(rawdata=data[,covs], data$treat, matched=matches$cnts, estimand='ATT')

    reg_form <- as.formula(data[, c('ppvtr.36', 'treat', covs)])
    reg_ps <- stan_glm(reg_form, data=data, weight=matches$cnts, algorithm='optimizing')

    print(summary(reg_ps)['treat', 1:2])

    return(list(pscores=pscores, matches=matches, bal=bal, reg_ps=reg_ps, covs=covs))
}


# Plot overlap of pscore
plot.overlap <- function(pscore, data=cc2){
    dev.off()
    par(mfrow=c(1,2))

    ctrl <- hist(pscore[data$treat==0])
    trt <- hist(pscore[data$treat==1])
    ymax <- round(max(c(ctrl$density, trt$density)), 1) + 1

    hist(pscore[data$treat==0], xlim=range(pscore), ylim=c(0,ymax),
        main='', border='darkgrey',
        mgp=c(2,.5,0), xlab='propensity scores', freq=FALSE)
    hist(pscore[data$treat==1], freq=FALSE, add=TRUE)

    ymax <- ceil(max(c(ctrl$counts, trt$counts))/10) * 10

    hist(pscore[data$treat==0], xlim=range(pscore), ylim=c(0,ymax),
        main='', border='darkgrey',
        mgp=c(2,.5,0), xlab='propensity scores', freq=TRUE)
    hist(pscore[data$treat==1], freq=TRUE, add=TRUE)
}


