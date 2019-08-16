# Generate code and figures for GHV CH 20

############################################

library(rstan)
    options(mc.cores=parallel::detectCores())
library(arm)
library(rstanarm)
library(survey)
source('library/matching.R')
source('library/balance.R')
source('library/estimation.R')

############################################

load('data/cc2.Rdata')


covs.all <- setdiff(names(cc2), c('row.names', 'row.names.1', 'treat', 'treat0', 'ppvtr.36'))
covs <- c('bw', 'preterm', 'dayskidh', 'sex', 'first', 'age', 'black', 'hispanic', 'white', 'b.marr', 'lths', 'hs', 'ltcoll', 'college', 'work.dur', 'prenatal', 'momage')
cov_names <- c('birth weight', 'weeks preterm', 'days in hospital', 'male', 'first born', 'age', 'black', 'hispanic', 'white', 'unmarried at birth', 'less than high school', 'high school graduate', 'some college', 'college graduate', 'worked during pregnancy', 'had no prenatal care', 'age at birth')

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


ps_spec <- formula(treat ~ bw + bwg + hispanic + black + b.marr + lths + hs + ltcoll + work.dur + prenatal + sex + first + preterm + momage + dayskidh + income)
ps_spec.st <- update(ps_spec, . ~ . + st5 + st9 + st12 + st25 + st36 + st42 + st48 + st53)

set.seed(20)
ps_fit_1 <- stan_glm(ps_spec, family=binomial(link='logit'), data=cc2, algorithm='optimizing')
ps_fit_1.st <- stan_glm(ps_spec.st, family=binomial(link='logit'), data=cc2, algorithm='optimizing')

pscores <- apply(posterior_linpred(ps_fit_1, type='link'), 2, mean)
pscores.st <- apply(posterior_linpred(ps_fit_1.st, type='link'), 2, mean)

# STEP 3: MATCHING
matches <- matching(z=cc2$treat, score=pscores, replace=FALSE)
matched <- cc2[matches$match.ind,]

matches.wr <- matching(z=cc2$treat, score=pscores, replace=TRUE)
wts.wr <- matches.wr$cnts

matches.st <- matching(z=cc2$treat, score=pscores.st, replace=FALSE)
matched.st <- cc2[matches.st$match.ind,]

matches.st.wr <- matching(z=cc2$treat, score=pscores.st, replace=TRUE)
wts.st.wr <- matches.st.wr$cnts

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
points(pts, idx, cex=1)
points(pts2, idx, pch=19, cex=1)
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
mar1 <- c(5, 4, 6, 2)
mar2 <- c(5, 3, 6, 4)
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
    main=main, cex.main=1.2,
    xlim=c(0,.55))
abline(v=0, lty=2)
points(pts, idx, cex=1)
points(pts2, idx, pch=19, cex=1)
axis(3, at=seq(0,.5,.1), xpd=TRUE)
axis(2, at=1:K, labels=longcovnames[1:K],
    las=2, hadj=1, lty=0, cex.axis=1)
# plot.balance(bal_nr.wr, longcovnames=cov_names, which.cov='cont', mar=c(1, 4, 5, 4))
pts <- bal_nr.wr$diff.means.raw[bal_nr.wr$binary==FALSE,4]
pts2 <- bal_nr.wr$diff.means.matched[bal_nr.wr$binary==FALSE,4]
# AZC: hack to fix spacing of binary covariates against x axis
# the spacing of how spaced apart the ticks are changes as the number of covariates change. It's frustratingly hard, maybe impossible, to get the spacing to match between the continuous and binary plots with different number of covariates in each, so, I'll add fake data that won't show up
pts <- c(pts, rep(NA, 7))
pts2 <- c(pts2, rep(NA, 7))
K <- length(pts)
idx <- 1:K
main <- 'Absolute Standardized Difference in Means'
longcovnames <- cov_names[bal_nr.wr$binary==FALSE]
# add extra names to match above
longcovnames <- c(longcovnames, rep('', 7))

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
segments(x0=0, y0=13, x1=0, y1=7.5, lty=2)
points(pts, idx, cex=1)
points(pts2, idx, pch=19, cex=1)
axis(3)
axis(2, at=8:12, labels=longcovnames[8:12],
    las=2, hadj=1, lty=0, cex.axis=1)
dev.off()
}


# Figure 20.14
############################################
# overlap of propensity scores before/after matching with replacement
{
pdf('outputs/ghv_ch20/ps.overlap.dens.AZC.pdf', width=11, height=8.5)
par(mfrow=c(1,2), cex.main=1.3, cex.lab=1.3)
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
par(mar=c(8,3,4,3), cex=1.4, cex.lab=1.2)
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
te_spec_nr <- update(ps_spec, ppvtr.36 ~ treat + .)

# treatment effect without replacement
set.seed(20)
reg_ps <- stan_glm(te_spec_nr, data=cc2[matches$match.ind,], algorithm='optimizing')
# treatment effect with replacement
set.seed(20)
reg_ps.wr <- stan_glm(te_spec_nr, data=cc2, weight=matches.wr$cnts, algorithm='optimizing')
ps_fit_1_design <- svydesign(ids=~1, weights=matches.wr$cnts, data=cc2)
reg_ps.wr_svy <- svyglm(te_spec_nr, design=ps_fit_1_design, data=cc2)

summary(reg_ps)['treat', 1:2]
summary(reg_ps.wr)['treat', 1:2]
summary(reg_ps.wr_svy)$coef['treat', 1:2]

# Geographic information, covs_nr.st
te_spec_nr.st <- update(ps_spec.st, ppvtr.36 ~ treat + . + st5 + st9 + st12 + st25 + st36 + st42 + st48 + st53)
# treatment effect without replacement
set.seed(20)
reg_ps.st <- stan_glm(te_spec_nr.st, data=cc2[matches.st$match.ind,], algorithm='optimizing')
# treatment effect with replacement
set.seed(20)
reg_ps.st.wr <- stan_glm(te_spec_nr.st, data=cc2, weight=matches.st.wr$cnts, algorithm='optimizing')
ps_fit_1.st_design <- svydesign(ids=~1, weights=matches.st.wr$cnts, data=cc2)
reg_ps.st.wr_svy <- svyglm(te_spec_nr.st, design=ps_fit_1.st_design, data=cc2)

summary(reg_ps.st)['treat', 1:2]
summary(reg_ps.st.wr)['treat', 1:2]
summary(reg_ps.st.wr_svy)$coef['treat', 1:2]



# standard regression estimate of treatment effect
set.seed(20)
reg_te <- stan_glm(te_spec_nr, data=cc2, algorithm='optimizing')
# standard regression estimate of treatment effect with state
set.seed(20)
reg_te.st <- stan_glm(te_spec_nr.st, data=cc2, algorithm='optimizing')

summary(reg_te)['treat', 1:2]
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

ps_spec2 <- formula(treat ~ bwg*as.factor(educ) + as.factor(ethnic)*b.marr + work.dur + prenatal + preterm + momage + sex + first + bw + dayskidT + pretermT + momageT + black*(bw + preterm + dayskidT) + b.marr*(bw + preterm + dayskidT) + bw*income)

set.seed(8)
ps_fit_2 <- stan_glm(ps_spec2, family=binomial(link="logit"), data=cc2, algorithm='optimizing')

pscores_2 <- apply(posterior_linpred(ps_fit_2, type='link'), 2, mean)
matches2 <- matching(z=cc2$treat, score=pscores_2, replace=FALSE)
matched2 <- cc2[matches2$match.ind,]
matches2_wr <- matching(z=cc2$treat, score=pscores_2, replace=TRUE)

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
plot.balance(bal_2.wr, longcovnames=cov_names)


# Treatment effect
te_spec2 <- formula(ppvtr.36 ~ treat + hispanic + black + b.marr + lths + hs + ltcoll + work.dur + prenatal + momage + sex + first + preterm + dayskidh + bw + income)
set.seed(8)
# MwoR
reg_ps2 <- stan_glm(te_spec2, data=matched2, algorithm='optimizing')
# MwR
reg_ps2.design <- svydesign(ids=~1, weights=matches2_wr$cnts, data=cc2)
reg_ps2.wr <- svyglm(te_spec2, design=reg_ps2.design, data=cc2)

summary(reg_ps2)['treat', 1:2]
summary(reg_ps2.wr)$coef['treat', 1:2]


# Geographic information using ps_spec2
ps_spec2.st <- update(ps_spec2, . ~ . + st5 + st9 + st12 + st25 + st36 + st42 + st48 + st53)

set.seed(8)
ps_fit_2.st <- stan_glm(ps_spec2.st, family=binomial(link='logit'), data=cc2, algorithm='optimizing')

pscores_2.st <- apply(posterior_linpred(ps_fit_2.st, type='link'), 2, mean)
matches2.st <- matching(z=cc2$treat, score=pscores_2.st, replace=FALSE)
matched2.st <- cc2[matches2.st$match.ind,]
matches2.st_wr <- matching(z=cc2$treat, score=pscores_2.st, replace=TRUE)

# Treatment effect estimate
te_spec2.st <- update(te_spec2, . ~ . + st5 + st9 + st12 + st25 + st36 + st42 + st48 + st53)
set.seed(8)
# MwoR
reg_ps2.st <- stan_glm(te_spec2.st, data=matched2.st, algorithm='optimizing')
reg_ps2.st_design <- svydesign(ids=~1, weights=matches2.st_wr$cnts, data=cc2)
reg_ps2.st.wr <- svyglm(te_spec2.st, design=reg_ps2.st_design, data=cc2)

summary(reg_ps2.st)['treat', 1:2]
summary(reg_ps2.st.wr)$coef['treat', 1:2]


############################################
# IPTW (including state indicators)

wt.iptw <- inv.logit(pscores) / (1 - inv.logit(pscores))
wt.iptw[cc2$treat==0] <- wt.iptw[cc2$treat==0] * (sum(wt.iptw[cc2$treat==0]) / sum(cc2$treat==0))
wt.iptw[cc2$treat==1] <- 1

set.seed(20)
ps_fit_iptw_design <- svydesign(ids=~1, weights=wt.iptw, data=cc2)
reg_ps.iptw <- svyglm(te_spec_nr, design=ps_fit_iptw_design, data=cc2)
summary(reg_ps.iptw)$coef['treat', 1:2]


############################################
# Section 20.8, Beyond balance in means
# table of ratio of standard deviations across treatment & control groups for unmatched, MWOR, MWR

cont_vars <- c('bw', 'preterm', 'dayskidh', 'age', 'momage', 'income')
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
# bw            0.50 0.91 0.52
# preterm       0.95 0.71 1.12
# dayskidh      2.07 0.82 3.13
# age           0.07 0.07 0.08
# momage        1.86 1.76 1.95
#          unmatched MWOR  MWR
# bw            0.50 0.90 0.52
# preterm       0.95 0.71 1.12
# dayskidh      2.07 0.82 2.76
# age           0.07 0.07 0.08
# momage        1.86 1.76 1.88
# income        0.27 0.19 0.29

# try with ps_fit_2
sds.mwor2 <- sapply(cont_vars, function(x){
    tapply(matched2[,x], matched2$treat, sd)
})
mwr.ind2 <- rep(cc2$row.names, times=matches2_wr$cnts)
matched2_wr <- cc2[mwr.ind2, ]
sds.mwr2 <- sapply(cont_vars, function(x){
    tapply(matched2_wr[,x], matched2_wr$treat, sd)
})

sd.ratios2 <- lapply(list(sds.um, sds.mwor2, sds.mwr2), function(mat){
    apply(mat, 2, function(col){
        col[2] / col[1]
    })
})

sd.table2 <- round(data.frame(
    unmatched=sd.ratios2[[1]],
    MWOR=sd.ratios2[[2]],
    MWR=sd.ratios2[[3]]), 2)
#          unmatched MWOR  MWR
# bw            0.50 0.84 0.59
# preterm       0.95 0.82 1.20
# dayskidh      2.07 0.88 3.07
# age           0.07 0.07 0.08
# momage        1.86 1.63 1.89
#          unmatched MWOR  MWR
# bw            0.50 0.85 0.56
# preterm       0.95 0.82 1.08
# dayskidh      2.07 0.88 2.82
# age           0.07 0.07 0.08
# momage        1.86 1.64 1.98
# income        0.27 0.23 0.78

# CBPS
library(CBPS)

psfit_cbps <- CBPS(treat ~ bwg*as.factor(educ) + as.factor(ethnic)*b.marr + work.dur + prenatal + preterm + momage + sex + first + bw + dayskidT +preterm + pretermT + momage + momageT + black*(bw + preterm + dayskidT) + b.marr*(bw + preterm + dayskidT) + bw*income, data=cc2, ATT=1, method='over')
pscores_cbps <- psfit_cbps$fitted.values
matches_cbps <- matching(cc2$treat, score=pscores_cbps, replace=FALSE)
matched_cbps <-  cc2[matches_cbps$match.ind,]
matches_wr_cbps <- matching(cc2$treat, score=pscores_cbps, replace=TRUE)
matched_wr_cbps <- cc2[matches_cbps$match.ind,]

sds.mwor_cbps <- sapply(cont_vars, function(x){
    tapply(matched_cbps[,x], matched_cbps$treat, sd)
})
sds.mwr_cbps <- sapply(cont_vars, function(x){
    tapply(matched_wr_cbps[,x], matched_wr_cbps$treat, sd)
})

sd.ratios3 <- lapply(list(sds.um, sds.mwor_cbps, sds.mwr_cbps), function(mat){
    apply(mat, 2, function(col){
        col[2] / col[1]
    })
})

sd.table3 <- round(data.frame(
    unmatched=sd.ratios3[[1]],
    MWOR=sd.ratios3[[2]],
    MWR=sd.ratios3[[3]]), 2)
#          unmatched MWOR  MWR
# bw            0.50 0.59 0.59
# preterm       0.95 0.75 0.75
# dayskidh      2.07 0.82 0.82
# age           0.07 0.06 0.06
# momage        1.86 1.29 1.29
#          unmatched MWOR  MWR
# bw            0.50 0.60 0.60
# preterm       0.95 0.74 0.74
# dayskidh      2.07 0.81 0.81
# age           0.07 0.06 0.06
# momage        1.86 1.27 1.27
# income        0.27 0.19 0.19

# genetic matching
mgen_1 <- readRDS('models/mgen_1.rds')
mgen_1.wr <- readRDS('models/mgen_1_wr.rds')
matched_mgen <- cc2[c(mgen_1$matches[,1], mgen_1$matches[,2]),]
matched_wr_mgen <- cc2[c(mgen_1.wr$matches[,1], mgen_1.wr$matches[,2]),]

sds.mwor_mgen <- sapply(cont_vars, function(x){
    tapply(matched_mgen[,x], matched_mgen$treat, sd)
})
sds.mwr_mgen <- sapply(cont_vars, function(x){
    tapply(matched_wr_mgen[,x], matched_wr_mgen$treat, sd)
})

sd.ratios4 <- lapply(list(sds.um, sds.mwor_mgen, sds.mwr_mgen), function(mat){
    apply(mat, 2, function(col){
        col[2] / col[1]
    })
})

sd.table4 <- round(data.frame(
    unmatched=sd.ratios4[[1]],
    MWOR=sd.ratios4[[2]],
    MWR=sd.ratios4[[3]]), 2)
#          unmatched MWOR  MWR
# bw            0.50 0.65 0.74
# preterm       0.95 0.77 0.90
# dayskidh      2.07 0.81 0.65
# age           0.07 0.08 0.09
# momage        1.86 1.70 1.88

mgen_2 <- readRDS('models/mgen_2.rds')
mgen_2.wr <- readRDS('models/mgen_2_wr.rds')
matched_mgen2 <- cc2[c(mgen_2$matches[,1], mgen_2$matches[,2]),]
matched_wr_mgen2 <- cc2[c(mgen_2.wr$matches[,1], mgen_2.wr$matches[,2]),]

sds.mwor_mgen2 <- sapply(cont_vars, function(x){
    tapply(matched_mgen2[,x], matched_mgen2$treat, sd)
})
sds.mwr_mgen2 <- sapply(cont_vars, function(x){
    tapply(matched_wr_mgen2[,x], matched_wr_mgen2$treat, sd)
})

sd.ratios4.2 <- lapply(list(sds.um, sds.mwor_mgen2, sds.mwr_mgen2), function(mat){
    apply(mat, 2, function(col){
        col[2] / col[1]
    })
})

sd.table4.2 <- round(data.frame(
    unmatched=sd.ratios4.2[[1]],
    MWOR=sd.ratios4.2[[2]],
    MWR=sd.ratios4.2[[3]]), 2)
