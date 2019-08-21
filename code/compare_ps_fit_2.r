############################################
# Compare ps_fits for GHV

############################################

library(rstan)
library(arm)
library(rstanarm)
    options(mc.cores=parallel::detectCores())
library(survey)
source('library/matching.R')
source('library/balance.R')
source('library/estimation.R')
source('code/library.r')

############################################

load('data/cc2.Rdata')

cc2$bwT = (cc2$bw-1500)^2
cc2$dayskidT = log(cc2$dayskidh+1)
cc2$pretermT = (cc2$preterm+8)^2
cc2$momageT = (cc2$momage^2)

covs.all <- setdiff(names(cc2), c('row.names', 'row.names.1', 'treat', 'treat0', 'ppvtr.36'))
covs <- c('bw', 'preterm', 'dayskidh', 'sex', 'first', 'age', 'black', 'hispanic', 'white', 'b.marr', 'lths', 'hs', 'ltcoll', 'college', 'work.dur', 'prenatal', 'momage')
cov_names <- c('birth weight', 'weeks preterm', 'days in hospital', 'male', 'first born', 'age', 'black', 'hispanic', 'white', 'unmarried at birth', 'less than high school', 'high school graduate', 'some college', 'college graduate', 'worked during pregnancy', 'had no prenatal care', 'age at birth')

############################################
# ps_fit_1
############################################
covs.nr <- c('bwg', 'hispanic', 'black', 'b.marr', 'lths', 'hs', 'ltcoll', 'work.dur', 'prenatal', 'sex', 'first', 'bw', 'preterm', 'momage', 'dayskidh')
covs.nr.st <- c(covs.nr, 'st5', 'st9', 'st12', 'st25', 'st36', 'st42', 'st53')

ps_spec <- formula(treat ~ bw + bwg + hispanic + black + b.marr + lths + hs + ltcoll + work.dur + prenatal + sex + first + preterm + momage + dayskidh + income)

set.seed(20)
ps_fit_1 <- stan_glm(ps_spec, family=binomial(link='logit'), data=cc2, algorithm='optimizing')

pscores <- apply(posterior_linpred(ps_fit_1, type='link'), 2, mean)

# STEP 3: MATCHING
matches <- matching(z=cc2$treat, score=pscores, replace=FALSE)
matched <- cc2[matches$match.ind,]

matches.wr <- matching(z=cc2$treat, score=pscores, replace=TRUE)
wts.wr <- matches.wr$cnts
matched.wr <- cc2[matches.wr$match.ind,]

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

############################################
# ps_fit_2
############################################

# Original ps_spec2 MwoR: 9.8 (1.5), MwR: 8.7 (2.4)
ps_spec2 <- formula(treat ~ bwg*as.factor(educ) + as.factor(ethnic)*b.marr + work.dur + prenatal + preterm + momage + sex + first + bw + dayskidT + pretermT + momageT + black*(bw + preterm + dayskidT) + b.marr*(bw + preterm + dayskidT) + bw*income)

set.seed(8)
ps_fit_2 <- stan_glm(ps_spec2, family=binomial(link="logit"), data=cc2, algorithm='optimizing')

pscores_2 <- apply(posterior_linpred(ps_fit_2, type='link'), 2, mean)
matches2 <- matching(z=cc2$treat, score=pscores_2, replace=FALSE)
matched2 <- cc2[matches2$match.ind,]
matches2_wr <- matching(z=cc2$treat, score=pscores_2, replace=TRUE)
matched2_wr <- cc2[matches2_wr$match.ind,]

te_spec2 <- formula(ppvtr.36 ~ treat + hispanic + black + b.marr + lths + hs + ltcoll + work.dur + prenatal + momage + sex + first + preterm + dayskidh + bw + income)
set.seed(8)
# MwoR
reg_ps2 <- stan_glm(te_spec2, data=matched2, algorithm='optimizing')
# MwR
reg_ps2.design <- svydesign(ids=~1, weights=matches2_wr$cnts, data=cc2)
reg_ps2.wr <- svyglm(te_spec2, design=reg_ps2.design, data=cc2)

summary(reg_ps2)['treat', 1:2]
summary(reg_ps2.wr)$coef['treat', 1:2]


############################################
# ps_fit_2 i21

# ps_specs_i21[242260] MwoR: 10.7 (1.5), MwR: 5.1 (2.6)
# ps_spec_i21 <- formula(treat~preterm+dayskidh+sex+first+black+hispanic+lths+hs+ltcoll+work.dur+prenatal+momage+income+dayskidT+pretermT+income+black:preterm+black:dayskidT+b.marr:bw+b.marr:preterm+b.marr:dayskidT)
# ps_specs_i21 [242260] MwoR: 10.7 (1.5), MwR: 5.1 (2.6)
# ps_spec_i21 <- formula(treat~preterm+dayskidh+sex+first+black+hispanic+lths+hs+ltcoll+work.dur+prenatal+momage+income+dayskidT+pretermT+income+black:preterm+black:dayskidT+b.marr:bw+b.marr:preterm+b.marr:dayskidT)
# ps_specs_i21[257460] MwoR: 11.0 (1.5), MwR: 8.1 (2.2)
# winner?
ps_spec_i21 <- formula(treat~preterm+dayskidh+sex+black+hispanic+b.marr+lths+hs+ltcoll+work.dur+prenatal+income+dayskidT+pretermT+momageT+income+black:preterm+black:dayskidT+b.marr:bw+b.marr:preterm+b.marr:dayskidT)
# ps_specs_i21 [95249] mwr balanced 20 covs. MwoR: 10.0 (1.5), MwR: 10.6 (2.3)
# ps_spec_i21 <- formula(treat~bw+preterm+dayskidh+sex+hispanic+b.marr+lths+hs+ltcoll+work.dur+prenatal+momage+income+bwT+pretermT+income+black:dayskidT+b.marr:bw+b.marr:preterm+b.marr:dayskidT+bw:income)

set.seed(8)
ps_fit_i21 <- stan_glm(ps_spec_i21, family=binomial(link="logit"), data=cc2, algorithm='optimizing')

pscores_i21 <- apply(posterior_linpred(ps_fit_i21, type='link'), 2, mean)
matches_i21 <- matching(z=cc2$treat, score=pscores_i21, replace=FALSE)
matched_i21 <- cc2[matches_i21$match.ind,]
matches_i21_wr <- matching(z=cc2$treat, score=pscores_i21, replace=TRUE)
matched_i21_wr <- cc2[matches_i21_wr$match.ind,]


set.seed(8)
# MwoR
reg_ps_i21 <- stan_glm(te_spec2, data=matched_i21, algorithm='optimizing')
# MwR
reg_ps_i21.design <- svydesign(ids=~1, weights=matches_i21_wr$cnts, data=cc2)
reg_ps_i21.wr <- svyglm(te_spec2, design=reg_ps_i21.design, data=cc2)

summary(reg_ps_i21)['treat', 1:2]
summary(reg_ps_i21.wr)$coef['treat', 1:2]

############################################
# ps_fit_2 i22

# ps_specs_i22[70824] winner? MwoR: 10.8 (1.6), MwR: 8.1 (2.9)
ps_spec_i22 <- formula(treat ~ preterm +dayskidh +sex +first +black +b.marr +lths +hs +ltcoll +work.dur +prenatal +momage +income +dayskidT +pretermT +income +black:preterm +black:dayskidT +b.marr:bw +b.marr:preterm +b.marr:dayskidT +bw:income)

set.seed(8)
ps_fit_i22 <- stan_glm(ps_spec_i22, family=binomial(link="logit"), data=cc2, algorithm='optimizing')

pscores_i22 <- apply(posterior_linpred(ps_fit_i22, type='link'), 2, mean)
matches_i22 <- matching(z=cc2$treat, score=pscores_i22, replace=FALSE)
matched_i22 <- cc2[matches_i22$match.ind,]
matches_i22_wr <- matching(z=cc2$treat, score=pscores_i22, replace=TRUE)
matched_i22_wr <- cc2[matches_i22_wr$match.ind,]

set.seed(8)
# MwoR
reg_ps_i22 <- stan_glm(te_spec2, data=matched_i22, algorithm='optimizing')
# MwR
reg_ps_i22.design <- svydesign(ids=~1, weights=matches_i22_wr$cnts, data=cc2)
reg_ps_i22.wr <- svyglm(te_spec2, design=reg_ps_i22.design, data=cc2)

summary(reg_ps_i22)['treat', 1:2]
summary(reg_ps_i22.wr)$coef['treat', 1:2]

############################################
# ps_fit_2 i23

# ps_specs_i23[1439] MwoR: 9.9 (1.5), MwR: 10. (2.8)
ps_spec_i23 <- formula(treat~bw+preterm+dayskidh+sex+first+black+hispanic+b.marr+lths+hs+ltcoll+prenatal+momage+income+bwT+dayskidT+pretermT+income+black:preterm+black:dayskidT+b.marr:bw+b.marr:dayskidT+bw:income)

set.seed(8)
ps_fit_i23 <- stan_glm(ps_spec_i23, family=binomial(link="logit"), data=cc2, algorithm='optimizing')

pscores_i23 <- apply(posterior_linpred(ps_fit_i23, type='link'), 2, mean)
matches_i23 <- matching(z=cc2$treat, score=pscores_i23, replace=FALSE)
matched_i23 <- cc2[matches_i23$match.ind,]
matches_i23_wr <- matching(z=cc2$treat, score=pscores_i23, replace=TRUE)
matched_i23_wr <- cc2[matches_i23_wr$match.ind,]

set.seed(8)
# MwoR
reg_ps_i23 <- stan_glm(te_spec2, data=matched_i23, algorithm='optimizing')
# MwR
reg_ps_i23.design <- svydesign(ids=~1, weights=matches_i23_wr$cnts, data=cc2)
reg_ps_i23.wr <- svyglm(te_spec2, design=reg_ps_i23.design, data=cc2)

summary(reg_ps_i23)['treat', 1:2]
summary(reg_ps_i23.wr)$coef['treat', 1:2]

############################################
# cobalt

bal_ihdp(matched)
bal_ihdp(matched2)
bal_ihdp(matched_i21)
bal_ihdp(matched_i22)
bal_ihdp(matched_i23)

bal_ihdp(matched.wr)
bal_ihdp(matched2_wr)
bal_ihdp(matched_i21_wr)
bal_ihdp(matched_i22_wr)
bal_ihdp(matched_i23_wr)
