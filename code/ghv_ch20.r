# Generate code and figures for GHV CH 20
# Only reproducing code and figures dependent on change in propensity score model, covariate selection, treatment effect estimation
############################################

library(dbarts)
library(rstan)
    options(mc.cores=parallel::detectCores())
library(rstanarm)
library(arm)
source('code/matching.R')
source('code/balance.R')
source('code/estimation.R')

############################################

load('data/cc2.Rdata')


# FIGURE 20.9
# balance plot
covs.all <- setdiff(names(cc2), c('row.names', 'row.names.1', 'treat', 'treat0', 'ppvtr.36'))
covs <- c('neg.bw', 'preterm', 'dayskidh', 'sex', 'first', 'age', 'black', 'hispanic', 'white', 'b.marr', 'lths', 'hs', 'ltcoll', 'college', 'work.dur', 'prenatal', 'momage')
cov_names <- c('negative birth weight', 'weeks preterm', 'days in hospital', 'male', 'first born', 'age', 'black', 'hispanic', 'white', 'unmarried at birth', 'less than high school', 'high school graduate', 'some college', 'college graduate', 'worked during pregnancy', 'had no prenatal care', 'age at birth')

form_20.9 <- as.formula(cc2[, c("treat", covs)])
ps_fit_20.9 <- stan_glm(form_20.9, family=binomial(link='logit'), data=cc2, algorithm='optimizing')

pscores_20.9 <- apply(posterior_linpred(ps_fit_20.9, type='link'), 2, mean)
matches_20.9 <- matching(z=cc2$treat, score=pscores_20.9, replace=FALSE)
matches_20.9.wr <- matching(z=cc2$treat, score=pscores_20.9, replace=TRUE)
bal_20.9 <- balance(rawdata=cc2[,covs], cc2$treat, matched=matches_20.9$cnts, estimand='ATT')
bal_20.9.wr <- balance(rawdat=cc2[,covs], cc2$treat, matched=matches_20.9.wr$cnts, estimand='ATT')

plot.balance(bal_20.9, longcovnames=cov_names)
plot.balance(bal_20.9.wr, longcovnames=cov_names)


# FIGURE 20.11
# overlap plot
pdf('outputs/ghv_ch20/age.educ.freq.AZC.f20.11.pdf', height=5, width=6)
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
# these are the no redundancy covariates with and without state covaraites
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

# STEP 4: DIAGNOSTICS FOR BALANCE AND OVERLAP
# separate balance plots for continuous and binary variables
pdf('outputs/ghv_ch20/balance.cont.binary.AZC.pdf', height=5, width=6)
par(mfrow=c(2,1))
plot.balance(bal_20.9.wr, longcovnames=cov_names, which.cov='cont')
plot.balance(bal_20.9.wr, longcovnames=cov_names, which.cov='binary')
dev.off()

# example: good overlap, bad pscore
ps3.mod <- glm(treat ~ unemp.rt, data=cc2,family=binomial) 
pscores3 <- predict(ps3.mod, type="link")

pdf('outputs/ghv_ch20/bad.pscore.overlap.AZC.pdf', width=4, height=4)
par(mfrow=c(1,2))
# Plot the overlapping histograms for pscore1b, density
hist(pscores3[cc2$treat==0], xlim=range(pscores3), ylim=c(0,8),
     main="", border="darkgrey", 
     mgp=c(2,.5,0), xlab="propensity scores",freq=FALSE)
hist(pscores3[cc2$treat==1], freq=FALSE, add=TRUE)
# Plot the overlapping histograms for pscore1b, frequency
hist(pscores3[cc2$treat==0], xlim=range(pscores3), ylim=c(0,1300),
     main="", border="darkgrey", 
     mgp=c(2,.5,0), xlab="propensity scores",freq=TRUE)
hist(pscores3[cc2$treat==1], freq=TRUE, add=TRUE)
dev.off()

# STEP 5: ESTIMATING A TREATMENT EFFECT USING THE RESTRUCTURED DATA
# treatment effect without replacement
reg_ps <- stan_glm(ppvtr.36 ~ treat + hispanic + black + b.marr + lths +hs + ltcoll + work.dur + prenatal + momage + sex + first + preterm + age + dayskidh + bw, data=cc2[matches$match.ind,], algorithm='optimizing')
summary(reg_ps)['treat', 1:2]

# treatment effect with replacement

# standard regression estimate of treatment effect
reg_te <- stan_glm(ppvtr.36 ~ treat + hispanic + black + b.marr + lths +hs + ltcoll + work.dur + prenatal + momage + sex + first + preterm + age + dayskidh + bw, data=cc2, algorithm='optimizing')
summary(reg_te)['treat', 1:2]
