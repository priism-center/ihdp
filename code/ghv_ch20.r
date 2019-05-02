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
source('code/pscore.r')

############################################

load('data/cc2.Rdata')

set.seed(8)


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
covs.nr <- c('bwg', 'hispanic', 'black', 'b.marr', 'lths', 'hs', 'ltcoll', 'work.dur', 'prenatal', 'sex', 'first', 'bw', 'preterm', 'momage', 'dayskidh')
covs.nr.st <- c(covs.nr, 'st5', 'st9', 'st12', 'st25', 'st36', 'st42', 'st53')
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

bal_nr <- balance(rawdata=cc2[,covs.nr], treat=cc2$treat, matched=matches$cnts, estimand='ATT')
bal_nr.st <- balance(rawdata=cc2[,covs.nr.st], treat=cc2$treat, matched=matches.st$cnts, estimand='ATT')

# STEP 4: DIAGNOSTICS FOR BALANCE AND OVERLAP
# separate balance plots for continuous and binary variables
pdf('outputs/ghv_ch20/balance.cont.binary.AZC.pdf', height=3, width=5)
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
reg_ps.wr <- stan_glm(ppvtr.36 ~ treat + hispanic + black + b.marr + lths +hs + ltcoll + work.dur + prenatal + momage + sex + first + preterm + age + dayskidh + bw, data=cc2, weight=matches.wr$cnts, algorithm='optimizing')
summary(reg_ps.wr)['treat', 1:2]

# treatment effect with replacement

# standard regression estimate of treatment effect
reg_te <- stan_glm(ppvtr.36 ~ treat + hispanic + black + b.marr + lths +hs + ltcoll + work.dur + prenatal + momage + sex + first + preterm + age + dayskidh + bw, data=cc2, algorithm='optimizing')
summary(reg_te)['treat', 1:2]



############################################
# FINDING A IMPROVED PSCORE MODEL
# find an improved pscore model using an interaction or squared term
############################################

# No redundancy covariate set, without states
all.nr <- pscoreAll(covs.nr, form=NULL, matching=FALSE, data=cc2)
all.nr.mwr <- pscoreAll(covs.nr, form=NULL, matching=TRUE, data=cc2)

all.nr.st <- pscoreAll(covs.nr, form=NULL, matching=FALSE, data=cc2)
all.nr.st.mwr <-  pscoreAll(covs.nr, form=NULL, matching=TRUE, data=cc2)

# No redundancy covariate set, without states, ethnic:b.marr
form1 <- as.formula(cc2[, c('treat', covs.nr)])
form1 <- update.formula(form1, ~ . + bwg:b.marr)
all.int1  <-  pscoreAll(covs=covs.nr, form=form1, matching=FALSE, data=cc2)


# STEPWISE MODEl SElECTION
# interactions
model_initial <- glm(treat ~ 1, family=binomial(link='logit'), data=cc2[, c('treat', covs.nr)])
upper_int <- as.formula(paste('treat ~ (', paste(covs.nr, collapse = '+'), ')^2'))

model_int <- stepAIC(model_initial, scope=list(lower=~1, upper=upper_int), direction='forward', trace=FALSE)

psfit_int <- stan_glm(model_int$formula, family=binomial(link='logit'), data=cc2[,c('treat',covs.nr)])
pscores_int <- apply(posterior_linpred(psfit_int, type='link'), 2, mean)
matches_int <- matching(z=cc2$treat, score=pscores_int, replace=FALSE)
bal_int <- balance(rawdata=cc2[,covs.nr], cc2$treat, matched=matches_int$cnts, estimand='ATT')
reg_int <- stan_glm(update.formula(model_int$formula, ppvtr.36 ~ treat + .), data=cc2[matches_int$match.ind,], algorithm='optimizing')

# quadratic
covs.nr.bin <- c('bwg', 'hispanic',  'black', 'b.marr', 'lths', 'hs', 'ltcoll', 'work.dur', 'prenatal', 'sex', 'first')
covs.nr.cont <- c('bw', 'preterm', 'momage', 'dayskidh')
upper_quad <- as.formula(paste('treat ~ ', paste(covs.nr, collapse = '+'), ' + poly(bw,2) + poly(preterm,2) + poly(momage,2) + poly(dayskidh,2)'))

model_quad <- stepAIC(model_initial, scope=list(lower=~1, upper=upper_quad), direction='forward', trace=FALSE)

psfit_quad <- stan_glm(model_quad$formula, family=binomial(link='logit'), data=cc2[,c('treat',covs.nr)])
pscores_quad <- apply(posterior_linpred(psfit_quad, type='link'), 2, mean)
matches_quad <- matching(z=cc2$treat, score=pscores_quad, replace=FALSE)
bal_quad <- balance(rawdata=cc2[,covs.nr], cc2$treat, matched=matches_quad$cnts, estimand='ATT')
reg_quad <- stan_glm(as.formula(cc2[,c('ppvtr.36', 'treat', covs.nr)]), data=cc2[matches_quad$matches.ind,], algorithm='optimizing')

# child age^2
upper <- as.formula(cc2[,c('treat',covs.nr,'age')])
upper <- update.formula(upper, ~ . + poly(age,2))

psfit_age <- stan_glm(upper, family='binomial', data=cc2, algorithm='optimizing')

pscores_age <- apply(posterior_linpred(psfit_age, type='link'), 2, mean)
matches_age <- matching(z=cc2$treat, score=pscores_age, replace=FALSE)
matches_age.wr <- matching(z=cc2$treat, score=pscores_age, replace=TRUE)
bal_age <- balance(rawdata=cc2[,c(covs.nr,'age')], cc2$treat, matched=matches_age$cnts, estimand='ATT')
bal_age.wr <- balance(rawdata=cc2[,c(covs.nr,'age')], cc2$treat, matched=matches_age.wr$cnts, estimand='ATT')

reg_age <- stan_glm(update.formula(upper, ppvtr.36 ~ treat + .), data=cc2[matches_age$match.ind,], algorithm='optimizing')
reg_age.wr <- stan_glm(update.formula(upper, ppvtr.36 ~ treat + .), data=cc2, weights=matches_age.wr$cnts, algorithm='optimizing')

# bw^2
form <- as.formula(cc2[, c('treat', covs.nr)])
form <- update.formula(form, ~ . + poly(bw,2))

psfit_bw <- stan_glm(form, family='binomial', data=cc2, algorithm='optimizing')

pscores_bw <- apply(posterior_linpred(psfit_bw, type='link'), 2, mean)
matches_bw <- matching(z=cc2$treat, score=pscores_bw, replace=FALSE)
matches_bw.wr <- matching(z=cc2$treat, score=pscores_bw, replace=TRUE)
bal_bw <- balance(rawdata=cc2[,c(covs.nr)], cc2$treat, matched=matches_bw$cnts, estimand='ATT')
bal_bw.wr <- balance(rawdata=cc2[,c(covs.nr)], cc2$treat, matched=matches_bw.wr$cnts, estimand='ATT')

reg_bw <- stan_glm(update.formula(form, ppvtr.36 ~ treat + .), data=cc2[matches_bw$match.ind,], algorithm='optimizing')
reg_bw.wr <- stan_glm(update.formula(form, ppvtr.36 ~ treat + .), data=cc2, weight=matches_bw.wr$cnts, algorithm='optimizing')


# preterm:momage
form <- as.formula(cc2[, c('treat', covs.nr)])
form <- update.formula(form, ~ . + preterm:momage)

psfit_preterm_momage <- stan_glm(form, family='binomial', data=cc2, algorithm='optimizing')

pscores_preterm_momage <- apply(posterior_linpred(psfit_preterm_momage, type='link'), 2, mean)
matches_preterm_momage <- matching(z=cc2$treat, score=pscores_preterm_momage, replace=FALSE)
matches_preterm_momage.wr <- matching(z=cc2$treat, score=pscores_preterm_momage, replace=TRUE)
bal_preterm_momage <- balance(rawdata=cc2[,c(covs.nr)], cc2$treat, matched=matches_preterm_momage$cnts, estimand='ATT')
bal_preterm_momage.wr <- balance(rawdata=cc2[,c(covs.nr)], cc2$treat, matched=matches_preterm_momage.wr$cnts, estimand='ATT')

reg_preterm_momage <- stan_glm(update.formula(form, ppvtr.36 ~ treat + .), data=cc2[matches_preterm_momage$match.ind,], algorithm='optimizing')
reg_preterm_momage.wr <- stan_glm(update.formula(form, ppvtr.36 ~ treat + .), data=cc2, weight=matches_preterm_momage.wr$cnts, algorithm='optimizing')

# ps model from GHV
covs_2 <- c('bwg', 'educ', 'ethnic', 'b.marr', 'work.dur', 'prenatal', 'preterm', 'age', 'momage', 'dayskidh', 'sex', 'first', 'bw', 'black', 'preterm')
ps_fit_2 <- stan_glm(treat ~ bwg*as.factor(educ) + as.factor(ethnic)*b.marr + work.dur + prenatal + preterm + age + momage + dayskidh + sex + first + bw + black*(bw + preterm) +b.marr*(bw + preterm), family=binomial(link="logit"), data=cc2, algorithm='optimizing')

pscores_2 <- apply(posterior_linpred(ps_fit_2, type='link'), 2, mean)
matches_2 <- matching(z=cc2$treat, score=pscores_2, replace=FALSE)
matches_2.wr <- matching(z=cc2$treat, score=pscores_2, replace=TRUE)
bal_2 <- balance(rawdata=cc2[,covs_2], cc2$treat, matched=matches_2$cnts, estimand='ATT')
bal_2.wr <- balance(rawdata=cc2[,covs_2], cc2$treat, matched=matches_2.wr$cnts, estimand='ATT')

reg_2 <- stan_glm(ppvtr.36 ~ treat + bwg*as.factor(educ) + as.factor(ethnic)*b.marr + work.dur + prenatal + preterm + age + momage + dayskidh + sex + first + bw + black*(bw + preterm) +b.marr*(bw + preterm), data=cc2[matches_2$match.ind,], algorithm='optimizing')
reg_2.wr <- stan_glm(ppvtr.36 ~ treat + bwg*as.factor(educ) + as.factor(ethnic)*b.marr + work.dur + prenatal + preterm + age + momage + dayskidh + sex + first + bw + black*(bw + preterm) +b.marr*(bw + preterm), data=cc2, weight=matches_2.wr$cnts, algorithm='optimizing')


