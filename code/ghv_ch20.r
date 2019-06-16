# Generate code and figures for GHV CH 20
# Only reproducing code and figures dependent on change in propensity score model, covariate selection, treatment effect estimation
############################################

library(dbarts)
library(rstan)
    options(mc.cores=parallel::detectCores())
library(rstanarm)
library(arm)
library(survey)
source('code/matching.R')
source('code/balance.R')
source('code/estimation.R')

############################################

load('data/cc2.Rdata')


covs.all <- setdiff(names(cc2), c('row.names', 'row.names.1', 'treat', 'treat0', 'ppvtr.36'))
covs <- c('neg.bw', 'preterm', 'dayskidh', 'sex', 'first', 'age', 'black', 'hispanic', 'white', 'b.marr', 'lths', 'hs', 'ltcoll', 'college', 'work.dur', 'prenatal', 'momage')
cov_names <- c('negative birth weight', 'weeks preterm', 'days in hospital', 'male', 'first born', 'age', 'black', 'hispanic', 'white', 'unmarried at birth', 'less than high school', 'high school graduate', 'some college', 'college graduate', 'worked during pregnancy', 'had no prenatal care', 'age at birth')

# Figure 20.9: See ps_fit_1 MwoR, prior to step 4

# FIGURE 20.11
# overlap plot
pdf('outputs/ghv_ch20/age.educ.freq.AZC.f20.11.pdf', height=4, width=6)
par(mfrow=c(1,2))
hist(cc2$educ[cc2$treat==0],xlim=c(0,5),main="",border="darkgrey",breaks=c(.5,1.5,2.5,3.5,4.5),mgp=c(2,.5,0),xlab="mother's education",freq=TRUE)
hist(cc2$educ[cc2$treat==1],xlim=c(0,5),xlab="education",breaks=c(.5,1.5,2.5,3.5,4.5),freq=TRUE,add=T)
hist(cc2$age[cc2$treat==0],xlim=c(0,110),main="",xlab="age of child (months)",border="darkgrey",breaks=seq(0,110,10),mgp=c(2,.5,0),
     freq=TRUE)
hist(cc2$age[cc2$treat==1],xlim=c(0,110),xlab="",breaks=seq(0,110,10),freq=TRUE, add=T)
dev.off()


# FIGURE 20.12
# table of stratified treatment effect estimate table
edu <- list(cc2$lths, cc2$hs, cc2$ltcoll, cc2$college)
# mean difference
y.trt <- sapply(edu, FUN=function(x) {
    tapply(cc2$ppvtr.36, list(cc2$treat, x), mean)[4]
})
y.ctrl <- sapply(edu, FUN=function(x) {
    tapply(cc2$ppvtr.36, list(cc2$treat, x), mean)[3]
})
te <- y.trt - y.ctrl
# sample sizes
n.trt <- sapply(edu, FUN=function(x) {
    tapply(rep(1, nrow(cc2)), list(cc2$treat, x), sum)[4]
})
n.ctrl <- sapply(edu, FUN=function(x) {
    sum(tapply(rep(1, nrow(cc2)), list(cc2$treat, x), sum)[3])
})
# std errors
var.trt <- sapply(edu, FUN=function(x) {
    tapply(cc2$ppvtr.36, list(cc2$treat, x), var)[4]
})
var.ctrl <- sapply(edu, FUN=function(x) {
    tapply(cc2$ppvtr.36, list(cc2$treat, x), var)[3]
})
se <- sqrt(var.trt/n.trt + var.ctrl/n.ctrl)
tes <- matrix(c(te, n.trt, n.ctrl, se), nrow=4)
rownames(tes) <- c('lths', 'hs', 'ltcoll', 'college')
colnames(tes) <- c('te', 'n.trt', 'n.ctrl', 'se')
round(tes, 1)


# EQUATION 20.4
# ATE
(9.3* 1358 + 4.1 * 1820 + 7.9* 837 + 4.6* 366) / (1358+1820+837+366)
# standard error
(1.5^2* 1358^2 + 1.9^2 * 1738^2 + 2.4^2* 789^2 + 2.3^2 * 366^2) / ((1358+1820+837+366)^2)


# EQUATION 20.5
# ATT
round((9.3* 126 + 4.1 * 82 + 7.9* 48 + 4.6* 34) / (126+82+48+34), 1)
# standard error
(1.5^2* 126^2 + 1.9^2 * 82^2 + 2.4^2* 48^2 + 2.3^2 * 34^2) / ((126+82+48+34)^2)


# STEP 2: ESTIMATING THE PROPENSITY SCORE
# these are the no redundancy covariates with and without state covariates
covs.nr <- c('bwg', 'hispanic', 'black', 'b.marr', 'lths', 'hs', 'ltcoll', 'work.dur', 'prenatal', 'sex', 'first', 'bw', 'preterm', 'momage', 'dayskidh')
covs.nr.st <- c(covs.nr, 'st5', 'st9', 'st12', 'st25', 'st36', 'st42', 'st53')


set.seed(20)
ps_fit_1 <- stan_glm(treat ~ bwg + hispanic + black + b.marr + lths + hs + ltcoll + work.dur + prenatal + sex + first + bw + preterm + momage + dayskidh, family=binomial(link='logit'), data=cc2, algorithm='optimizing')
ps_fit_1.st <- stan_glm(treat ~ bwg + hispanic + black + b.marr + lths + hs + ltcoll + work.dur + prenatal + sex + first + st5 + st9 + st12 + st25 + st36 + st42 + st48 + st53 + bw + preterm + momage + dayskidh, family=binomial(link='logit'), data=cc2, algorithm='optimizing')

pscores <- apply(posterior_linpred(ps_fit_1, type='link'), 2, mean)
pscores.st <- apply(posterior_linpred(ps_fit_1.st, type='link'), 2, mean)

# STEP 3: MATCHING
matches <- matching(z=cc2$treat, score=pscores, replace=FALSE)
matched <- cc2[matches$match.ind,]

matches.wr <- matching(z=cc2$treat, score=pscores, replace=TRUE)
wts.wr <- matches.wr$cnts

matches.st <- matching(z=cc2$treat, score=pscores.st, replace=FALSE)
matched.st <- cc2[matches$match.ind,]

matches.wr.st <- matching(z=cc2$treat, score=pscores.st, replace=TRUE)
wts.wr.st <- matches.wr$cnts

# Balance plots for all covariattes, not just those used for pscore spec
bal_nr <- balance(rawdata=cc2[,covs], treat=cc2$treat, matched=matches$cnts, estimand='ATT')
bal_nr.wr <- balance(rawdata=cc2[,covs], treat=cc2$treat, matched=matches.wr$cnts, estimand='ATT')
bal_nr.st <- balance(rawdata=cc2[, union(covs, covs.nr.st)], treat=cc2$treat, matched=matches.st$cnts, estimand='ATT')

# Figure 20.9
############################################
# balance plot, labelled cov names, ps_fit_1 MwoR
# manual plot code taken from balance.R
{
pdf('outputs/ghv_ch20/balance.both.azc.pdf', width=11, height=8.5)
# balance for all covariates, not just those used in propensity score model
pts <- bal_nr$diff.means.raw[,4]
pts2 <- bal_nr$diff.means.matched[,4]
K <- length(pts)
idx <- 1:K
main <- 'Absolute Standardized Difference in Means'

mar <- c(18, 6, 6, 7)
par(mar=mar)

maxchar <- max(sapply(cov_names, nchar))
min.mar <- par('mar')
mar[2] <- max(min.mar[2], trunc(mar[2] + maxchar/10)) + mar[2] + 0.5
par(mar=mar)

pts <- rev(pts)
pts2 <- rev(pts2)
longcovnames <- rev(cov_names)

plot(c(pts,pts2), c(idx,idx),
    bty='n', xlab='', ylab='',
    xaxt='n', yaxt='n', type='n',
    main=main, cex.main=1.2)
abline(v=0, lty=2)
points(pts, idx, cex=1.8)
points(pts2, idx, pch=19, cex=1.8)
axis(3)
axis(2, at=1:K, labels=longcovnames[1:K],
    las=2, hadj=1, lty=0, cex.axis=1.2)
dev.off()
}
############################################

# STEP 4: DIAGNOSTICS FOR BALANCE AND OVERLAP
# separate balance plots for continuous and binary variables
# Figure 20.13

{
pdf('outputs/ghv_ch20/balance.cont.binary.AZC.pdf', width=10, height=6)
par(mfrow=c(1,2))
mar1 <- c(6, 4, 6, 4)
mar2 <- c(18, 3, 6, 4)
# plot.balance(bal_nr.wr, longcovnames=cov_names, which.cov='binary', mar=c(2, 4, 5, 4))
pts <- bal_nr.wr$diff.means.raw[bal_nr.wr$binary==TRUE,4]
pts2 <- bal_nr.wr$diff.means.matched[bal_nr.wr$binary==TRUE,4]
K <- length(pts)
idx <- 1:K
main <- 'Absolute Difference in Means'
longcovnames <- cov_names[bal_nr.wr$binary==TRUE]

mar <- mar1
par(mar=mar)
maxchar <- max(sapply(longcovnames, nchar))
min.mar <- par('mar')
mar[2] <- max(min.mar[2], trunc(mar[2] + maxchar/10)) + mar[2] + 0.5
par(mar=mar)

pts <- rev(pts)
pts2 <- rev(pts2)
longcovnames <- rev(longcovnames)

plot(c(pts,pts2), c(idx,idx),
    bty='n', xlab='', ylab='',
    xaxt='n', yaxt='n', type='n',
    main=main, cex.main=1.2)
abline(v=0, lty=2)
points(pts, idx, cex=1.6)
points(pts2, idx, pch=19, cex=1.6)
axis(3)
axis(2, at=1:K, labels=longcovnames[1:K],
    las=2, hadj=1, lty=0, cex.axis=1)
# plot.balance(bal_nr.wr, longcovnames=cov_names, which.cov='cont', mar=c(1, 4, 5, 4))
pts <- bal_nr.wr$diff.means.raw[bal_nr.wr$binary==FALSE,4]
pts2 <- bal_nr.wr$diff.means.matched[bal_nr.wr$binary==FALSE,4]
K <- length(pts)
idx <- 1:K
main <- 'Absolute Standardized Difference in Means'
longcovnames <- cov_names[bal_nr.wr$binary==FALSE]

mar <- mar2
par(mar=mar)
maxchar <- max(sapply(longcovnames, nchar))
min.mar <- par('mar')
mar[2] <- max(min.mar[2], trunc(mar[2] + maxchar/10)) + mar[2] + 0.5
par(mar=mar)

pts <- rev(pts)
pts2 <- rev(pts2)
longcovnames <- rev(longcovnames)

plot(c(pts,pts2), c(idx,idx),
    bty='n', xlab='', ylab='',
    xaxt='n', yaxt='n', type='n',
    main=main, cex.main=1.2)
abline(v=0, lty=2)
points(pts, idx, cex=1.6)
points(pts2, idx, pch=19, cex=1.6)
axis(3)
axis(2, at=1:K, labels=longcovnames[1:K],
    las=2, hadj=1, lty=0, cex.axis=1)
dev.off()
}


# Figure 20.14
############################################
# overlap of propensity scores before/after matching with replacement
{
pdf('outputs/ghv_ch20/ps.overlap.dens.AZC.pdf', width=11, height=8.5)
par(mfrow=c(1,2))
# Plot the overlapping histograms for pscores before matching, density
par(mar=c(16,8,2,2))
hist(pscores[cc2$treat==0], xlim=c(-20,5), ylim=c(0,.28), main="before matching", border="darkgrey", mgp=c(2,.5,0), xlab="logit propensity scores", freq=FALSE)
hist(pscores[cc2$treat==1], freq=FALSE, add=TRUE)
# Plot the overlapping histograms for pscores after matching, frequency
par(mar=c(16,3,2,8))
hist(pscores[cc2[matches.wr$match.ind, 'treat']==0], xlim=c(-20,6), ylim=c(0,.28), main="after matching", border="darkgrey", mgp=c(2,.5,0), xlab="logit propensity scores", freq=FALSE)
hist(pscores[cc2[matches.wr$match.ind, 'treat']==1], freq=FALSE, add=TRUE)
dev.off()
}

# how many pscores[cc2$treat==0] left out of plot?
sum(pscores[cc2$treat==0] < -20)

# pscore matching check
sum(pscores[cc2$treat==1] > max(pscores[cc2$treat==0]))


# Figures 20.16
############################################
# example: good overlap, bad pscore
set.seed(20)
ps3.mod <- glm(treat ~ unemp.rt, data=cc2,family=binomial) 
pscores3 <- predict(ps3.mod, type="link")

{
pdf('outputs/ghv_ch20/bad.pscore.overlap.AZC.pdf', width=11, height=8.5)
par(mar=c(8,3,4,3), cex=1.4)
# par(mar=c(16,8,2,2))
# Plot the overlapping histograms for pscore3, density
hist(pscores3[cc2$treat==0], xlim=range(pscores3), ylim=c(0,8),
     main="", border="darkgrey", 
     mgp=c(2,.5,0), xlab="logit propensity scores",freq=FALSE)
hist(pscores3[cc2$treat==1], freq=FALSE, add=TRUE)
dev.off()
}

############################################
# STEP 5: ESTIMATING A TREATMENT EFFECT USING THE RESTRUCTURED DATA
# treatment effect without replacement
set.seed(20)
reg_ps <- stan_glm(ppvtr.36 ~ treat + hispanic + black + b.marr + lths +hs + ltcoll + work.dur + prenatal + momage + sex + first + preterm + age + dayskidh + bw, data=cc2[matches$match.ind,], algorithm='optimizing')
summary(reg_ps)['treat', 1:2]
# treatment effect with replacement
reg_ps.wr <- stan_glm(ppvtr.36 ~ treat + hispanic + black + b.marr + lths +hs + ltcoll + work.dur + prenatal + momage + sex + first + preterm + age + dayskidh + bw, data=cc2, weight=matches.wr$cnts, algorithm='optimizing')
summary(reg_ps.wr)['treat', 1:2]
ps_fit_1_design <- svydesign(ids=~1, weights=matches.wr$cnts, data=cc2)
reg_ps.wr_svy <- svyglm(ppvtr.36 ~ treat + hispanic + black + b.marr + lths +hs + ltcoll + work.dur + prenatal + momage + sex + first + preterm + age + dayskidh + bw, design=ps_fit_1_design, data=cc2)
summary(reg_ps.wr_svy)$coef['treat', 1:2]

# covs_nr.st; state variables
set.seed(20)
reg_ps.st <- stan_glm(ppvtr.36 ~ treat + hispanic + black + b.marr + lths +hs + ltcoll + work.dur + prenatal + momage + sex + first + preterm + age + dayskidh + bw + st5 + st9 + st12 + st25 + st36 + st42 + st48 + st53, data=cc2[matches.st$match.ind,], algorithm='optimizing')
summary(reg_ps.st)['treat', 1:2]
# treatment effect with replacement
reg_ps.wr.st <- stan_glm(ppvtr.36 ~ treat + hispanic + black + b.marr + lths +hs + ltcoll + work.dur + prenatal + momage + sex + first + preterm + age + dayskidh + bw + st5 + st9 + st12 + st25 + st36 + st42 + st48 + st53, data=cc2, weight=matches.wr.st$cnts, algorithm='optimizing')
summary(reg_ps.wr.st)['treat', 1:2]
ps_fit_1_design.st <- svydesign(ids=~1, weights=matches.wr.st$cnts, data=cc2)
reg_ps.wr_svy.st <- svyglm(ppvtr.36 ~ treat + hispanic + black + b.marr + lths +hs + ltcoll + work.dur + prenatal + momage + sex + first + preterm + age + dayskidh + bw + st5 + st9 + st12 + st25 + st36 + st42 + st48 + st53, design=ps_fit_1_design.st, data=cc2)
summary(reg_ps.wr_svy.st)$coef['treat', 1:2]



# standard regression estimate of treatment effect
set.seed(20)
reg_te <- stan_glm(ppvtr.36 ~ treat + hispanic + black + b.marr + lths +hs + ltcoll + work.dur + prenatal + momage + sex + first + preterm + age + dayskidh + bw, data=cc2, algorithm='optimizing')
summary(reg_te)['treat', 1:2]

# standard regression estimate of treatment effect with state
reg_te.st <- stan_glm(ppvtr.36 ~ treat + hispanic + black + b.marr + lths +hs + ltcoll + work.dur + prenatal + momage + sex + first + preterm + age + dayskidh + bw + st5 + st9 + st12 + st25 + st36 + st42 + st48 + st53, data=cc2, algorithm='optimizing')
summary(reg_te.st)['treat', 1:2]



############################################
# FINDING A IMPROVED PSCORE MODEL
# find an improved pscore model using an interaction or squared term
############################################

# transformed variables
cc2$bwT = (cc2$bw-1500)^2
cc2$dayskidT = log(cc2$dayskidh+1)
cc2$pretermT = (cc2$preterm+8)^2
cc2$momageT = (cc2$momage^2)

set.seed(8)
ps_fit_2 <- stan_glm(treat ~ bwg*as.factor(educ) + as.factor(ethnic)*b.marr + work.dur + prenatal + preterm + age + momage + sex + first + bw + dayskidT +preterm + pretermT + momage + momageT + black*(bw + preterm + dayskidT) + b.marr*(bw + preterm + dayskidT), family=binomial(link="logit"), data=cc2, algorithm='optimizing')

pscores_2 <- apply(posterior_linpred(ps_fit_2, type='link'), 2, mean)
matches2 <- matching(z=cc2$treat, score=pscores_2, replace=FALSE)
matches2_wr <- matching(z=cc2$treat, score=pscores_2, replace=TRUE)
matched2_wr <- cc2[matches2_wr$matched,]

bal_2 <- balance(rawdata=cc2[,covs], cc2$treat, matched=matches2$cnts, estimand='ATT')
bal_2.wr <- balance(rawdata=cc2[,covs], cc2$treat, matched=matches2_wr$cnts, estimand='ATT')

par(mfrow=c(2,1))
plot.balance(bal_nr, which.covs='cont', main='ps_fit_1')
plot.balance(bal_2.wr, which.covs='cont', main='ps_fit_2')
plot.balance(bal_nr, which.covs='binary', main='ps_fit_1')
plot.balance(bal_2.wr, which.covs='binary', main='ps_fit_2')
dev.off()

# new figure
# side by side binary/continuous, ps_fit_2.wr
pdf
plot.balance(bal_2.wr, longcovnames=cov_names)

reg_ps <- stan_glm(ppvtr.36 ~ treat + hispanic + black + b.marr + lths + hs + ltcoll + work.dur + prenatal + momage + sex + first + preterm + age + dayskidh + bw, data=cc2, algorithm='optimizing')
summary(reg_ps)['treat', 1:2]

reg_ps2_design <- svydesign(ids=~1, weights=~matches2_wr$cnts, data=cc2)
reg_ps2.wr <- svyglm(ppvtr.36 ~ treat + hispanic + black + b.marr + lths + hs + ltcoll + work.dur + prenatal + momage + sex + first + preterm + age + dayskidh + bw, design=reg_ps2_design, data=cc2)
summary(reg_ps2)$coef['treat', 1:2]


############################################
# IPTW (including state indicators)

wt.iptw <- inv.logit(pscores) / (1 - inv.logit(pscores))
wt.iptw[cc2$treat==0] <- wt.iptw[cc2$treat==0] * (sum(wt.iptw[cc2$treat==0]) / sum(cc2$treat==0))
wt.iptw[cc2$treat==1] <- 1

ps_fit_iptw_design <- svydesign(ids=~1, weights=wt.iptw, data=cc2)
reg_ps.iptw <- svyglm(ppvtr.36 ~ treat + hispanic + black + b.marr + lths +hs + ltcoll + work.dur + prenatal + momage + sex + first + preterm + age + dayskidh + bw, design=ps_fit_iptw_design, data=cc2)
summary(reg_ps.iptw)$coef['treat', 1:2]


############################################
# Section 20.8, Beyond balance in means
# table of ratio of standard deviations across treatment & control groups for unmatched, MWOR, MWR

cont_vars <- c('neg.bw', 'preterm', 'dayskidh', 'age', 'momage')
sds.um <- sapply(cont_vars, function(x){
    tapply(cc2[,x], cc2$treat, sd)
})
sds.mwor <- sapply(cont_vars, function(x){
    tapply(matched[,x], matched$treat, sd)
})
mwr.ind <- rep(cc2$row.names, times=matches.wr$cnts)
matched.wr <- cc2[mwr.ind, ]
sds.mwr <- sapply(cont_vars, function(x){
    tapply(matched.wr[,x], matched.wr$treat, sd)
})

sd.ratios <- lapply(list(sds.um, sds.mwor, sds.mwr), function(mat){
    apply(mat, 2, function(col){
        col[2] / col[1]
    })
})

sd.table <- round(data.frame(
    unmatched=sd.ratios[[1]],
    MWOR=sd.ratios[[2]],
    MWR=sd.ratios[[3]]), 2)
#          unmatched MWOR  MWR
# neg.bw        0.50 0.91 0.52
# preterm       0.95 0.71 1.12
# dayskidh      2.07 0.82 3.13
# age           0.07 0.07 0.08
# momage        1.86 1.76 1.95

