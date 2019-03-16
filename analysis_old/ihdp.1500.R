### "truth" for this group is 7.4

# code for looking at nlsy/ihdp high quality child care data

cc2<-cc
cc2<-cc2[cc$bw>1500,]

cc2$neg.bw = 2500 - cc2$bw
cc2$no.prenatal = 1-cc2$prenatal
cc2$b.unmarr = 1-cc2$b.marr

covs <- c("hispanic","black","white","b.marr","lths","hs","ltcoll","college","work.dur","prenatal","momage","sex","first","preterm","age","dayskidh","bw")
covs.st <- c("hispanic","black","white","b.marr","lths","hs","ltcoll","college","work.dur","prenatal","momage","sex","first","age","preterm","dayskidh","bw","unemp.rt","st5","st9","st12","st25","st36","st42","st48","st53")

covs.nba <- c("hispanic","black","white","b.marr","lths","hs","ltcoll","college",
             "work.dur","prenatal","momage","sex","first","age")
covs.ba <- c("dayskidh","bw","preterm")

cc2$bwT = (cc2$bw-1500)^2
#cc2$bwT = (cc2$neg.bw+5100)^2
cc2$dayskidT = log(cc2$dayskidh+1)
cc2$pretermT = (cc2$preterm+8)^2
cc2$momageT = (cc2$momage^2)

save.image()
#### at this point save this dataset to a file
save(cc2,file="cc2.Rdata")

#######
## make a dataset to investigate bias amplification
data.ba = cc2[,c("ppvtr.36", "treat", covs.ba,covs.nba)]
lm(data.ba)

data.nba = cc2[,c("ppvtr.36", "treat",covs.nba)]
lm(data.nba)

library(foreign)
write.dta(data.ba,file="data_ba.dta")
#######

psdata = cc2[,covs.st]
save(psdata,file="psdata.Rdata")

#########  FIRST OVERLAP PLOTS


#### UNIVARIATE
postscript("bw.overlap.dens.ps", height=3.5in, width=5in, horizontal=T)
par(mfrow=c(1,1))
hist(cc2$bw[cc2$treat==1],xlim=c(0,6000),main="treatment group",xlab="birthweight",freq=F)
hist(cc2$bw[cc2$treat==0],xlim=c(0,6000),main="control group",xlab="birthweight",freq=F,add=T,border="darkgrey")
dev.off()

postscript("bw.overlap.freq.ps", height=3.5in, width=5in, horizontal=T)
par(mfrow=c(1,1))
hist(cc2$bw[cc2$treat==0],xlim=c(0,6000),main="control group",xlab="birthweight",border="darkgrey")
hist(cc2$bw[cc2$treat==1],xlim=c(0,6000),main="treatment group",xlab="birthweight",add=T)
dev.off()

# breaks the same across groups
postscript("bw.overlap.freq.ps", height=3.5in, width=5in, horizontal=T)
par(mfrow=c(1,1))
hist(cc2$bw[cc2$treat==0],xlim=c(0,6000),main="control group",xlab="birthweight",border="darkgrey",breaks=seq(0,6000,400))
hist(cc2$bw[cc2$treat==1],xlim=c(0,6000),main="treatment group",xlab="birthweight",add=T,breaks=seq(0,6000,400))
dev.off()


postscript("age.overlap.dens.ps", height=3.5in, width=5in, horizontal=T)
par(mfrow=c(1,1))
hist(cc2$age[cc2$treat==1],xlim=c(0,110),main="treatment group",xlab="age",breaks=seq(0,110,10),freq=F)
hist(cc2$age[cc2$treat==0],xlim=c(0,110),main="control group",xlab="age",border="darkgrey",breaks=seq(0,110,10),add=T,freq=F)
dev.off()


## now preterm
postscript("preterm.overlap.freq.ps", height=3.5in, width=5in, horizontal=T)
par(mfrow=c(1,1))
hist(cc2$preterm[cc2$treat==0],xlim=c(0,6000),main="control group",xlab="birthweight",border="darkgrey")
hist(cc2$preterm[cc2$treat==1],xlim=c(0,6000),main="treatment group",xlab="birthweight",add=T)
dev.off()

## now educ

par(mfrow=c(1,1))
hist(cc$educ[cc$treat==0],xlim=c(0,5),main="control group",xlab="education",border="darkgrey",freq=F)
hist(cc$educ[cc$treat==1],xlim=c(0,5),main="treatment group",xlab="education",add=T,freq=F)

par(mfrow=c(1,1))
hist(cc$college[cc$treat==0],xlim=c(0,5),main="control group",xlab="college",border="darkgrey",freq=F)
hist(cc$college[cc$treat==1],xlim=c(0,5),main="treatment group",xlab="college",add=T,freq=F)

par(mfrow=c(1,1))
hist(cc$educ3[cc$treat==0],xlim=c(0,5),main="control group",xlab="educ3",border="darkgrey",freq=F)
hist(cc$educ3[cc$treat==1],xlim=c(0,5),main="treatment group",xlab="educ3",add=T,freq=F)

## plots actually used to illustrate distinctions between overlap and balance
# now age
postscript("age.educ.freq.ps", height=4, width=6, horizontal=T)
par(mfrow=c(1,2))
plot(x=seq(0,5,.1),y=seq(0,1800,(1800/50)),bty="n",xaxt="n",yaxt="n",mgp=c(2,.5,0),xlab="mother's education",ylab="frequency",type="n",main="")
axis(1, 1:4)
axis(2, c(0,1000,2000))
hist(cc2$educ[cc$treat==0],xlim=c(0,5),main="",border="darkgrey",breaks=c(.5,1.5,2.5,3.5,4.5),add=T)
hist(cc2$educ[cc$treat==1],xlim=c(0,5),xlab="education",breaks=c(.5,1.5,2.5,3.5,4.5),add=T)
#
hist(cc2$age[cc2$treat==0],xlim=c(0,110),main="",xlab="age of child (months)",border="darkgrey",breaks=seq(0,110,10),mgp=c(2,.5,0))
hist(cc2$age[cc2$treat==1],xlim=c(0,110),xlab="",breaks=seq(0,110,10),add=T)
dev.off()

#covs <- c("hispanic","black","white","b.marr","lths","hs","ltcoll","college","work.dur","prenatal","momage","sex","first","preterm","age","dayskidh","bw","unemp.rt","st5","st9","st12","st25","st36","st42","st48","st53")
covs2 <- c("neg.bw","preterm","dayskidh","sex","first","age","black","hispanic","white","b.unmarr","lths","hs","ltcoll","college","work.dur","no.prenatal","momage") 
cov.nms <- c("negative birth weight","weeks preterm","days in hospital","male","first born","age","black","Hispanic","white","unmarried at birth","less than high school","high school graduate","some college","college graduate","worked during pregnancy","had no prenatal care","age at birth")

diff.means=matrix(0,length(covs2),6)
for(i in 1:length(covs2)){
diff.means[i,1:2] <- c(mean(cc2[cc2$treat==1,covs2[i]]),mean(cc2[cc2$treat==0,covs2[i]]))
diff.means[i,3] <- diff.means[i,1]-diff.means[i,2]
diff.means[i,5] <- sqrt(var(cc2[cc2$treat==1,covs2[i]])/sum(cc2$treat==1) + var(cc2[cc2$treat==0,covs2[i]])/sum(cc2$treat==0))
diff.means[i,6] <- sqrt((var(cc2[cc2$treat==1,covs2[i]])+var(cc2[cc2$treat==0,covs2[i]]))/2)
diff.means[i,4] <- diff.means[i,3]/diff.means[i,6]
}
dimnames(diff.means) <- list(covs2,c("treat","control","diff","diff.std","se","sd"))
round(diff.means,2)

regression.simp.table(name=cov.nms, est=diff.means[,3], sd=diff.means[,6], file="balance.table.ps")

           treat control     diff diff.std    se     sd
neg.bw      491.35 -835.27 1326.62     2.98 18.81 444.70
preterm       6.07    1.18    4.89     2.49  0.12   1.97
dayskidh     14.69    4.17   10.52     1.19  0.67   8.86
sex           0.51    0.50    0.01     0.02  0.03   0.50
first         0.48    0.42    0.06     0.12  0.03   0.50
age          56.67   56.35    0.32     0.02  0.47  20.59
black         0.50    0.28    0.22     0.46  0.03   0.48
hispanic      0.09    0.21   -0.12    -0.34  0.02   0.36
white         0.40    0.50   -0.10    -0.20  0.03   0.50
b.unmarr      0.57    0.31    0.26     0.53  0.03   0.48
lths          0.43    0.30    0.13     0.28  0.03   0.48
hs            0.28    0.42   -0.14    -0.30  0.03   0.47
ltcoll        0.17    0.19   -0.03    -0.07  0.02   0.38
college       0.12    0.08    0.04     0.12  0.02   0.30
work.dur      0.59    0.62   -0.03    -0.06  0.03   0.49
no.prenatal   0.04    0.01    0.03     0.19  0.01   0.17
momage       24.44   23.75    0.69     0.15  0.35   4.71


#### BIVARIATE

#need to jitter things a little

# birthweight
postscript("ppvt.bw.ps", horizontal=T,height=3.6,width=4.2)
par(mfrow=c(1,1))
tmp <- lm(ppvtr.36~bw+treat,data=cc2)$coef
#plot(cc2$bw, cc2$ppvtr.36, xlab="birth weight", ylab="test score at age 3", cex=.8, mgp=c(2,.5,0), main="",type="n",xlim=c(1500,5000))
plot(cc2$bw, cc2$ppvtr.36, xlab="birth weight", ylab="test score at age 3", mgp=c(2,.5,0), main="",type="n",xlim=c(1500,5000),cex.axis=.75,cex.lab=.8,lab=c(3,5,7),xaxt="n")
axis(side=1,at=c(2000,3000,4000,5000),cex.axis=.75)
points(cc2$bw[cc2$treat==0]+runif(sum(cc2$treat==0),-.5,5), cc2$ppvtr.36[cc2$treat==0], col="darkgrey", pch=20, cex=.3)
points(cc2$bw[cc2$treat==1]+runif(sum(cc2$treat==1),-.5,5), cc2$ppvtr.36[cc2$treat==1], pch=20, cex=.3)
#curve(tmp[1]+tmp[2]*x,add=T,col="darkgrey",cex=2)
curve(tmp[1]+tmp[2]*x,add=T,lty=2)
curve(tmp[1]+tmp[3]+tmp[2]*x,add=T)
dev.off()
## also tried ineeds, age, momage (shows treats with no comps), momed (excellent
##overlap)

# income to needs
postscript("ppvt.age.ps", horizontal=T)
par(mfrow=c(1,1))
tmp <- lm(ppvtr.36~age+treat,data=cc22)$coef
plot(cc2$age, cc2$ppvtr.36, xlab="age of child (months) in year 2000", ylab="test scores at age 3", cex=.8, mgp=c(2,.5,0), main="data and estimated regression line",type="n")
points(cc2$age[cc2$treat==1], cc2$ppvtr.36[cc2$treat==1], pch=20, cex=.8)
points(cc2$age[cc2$treat==0], cc2$ppvtr.36[cc2$treat==0], cex=.8, pch=20, col="gray")
curve(tmp[1]+tmp[2]*x,add=T)
curve(tmp[1]+tmp[3]+tmp[2]*x,add=T)
dev.off()

####### NOW STRATIFICATION

tes.ns=matrix(0,4,4)
tes.ns[1,]<-c(mean(cc2$ppvtr.36[cc2$treat==1 & cc2$lths==1]) - mean(cc2$ppvtr.36[cc2$treat==0 & cc2$lths==1]),sum(cc2$treat==1 & cc2$lths==1),sum(cc2$lths==1),var(cc2$ppvtr.36[cc2$treat==1 & cc2$lths==1])/sum(cc2$treat==1 & cc2$lths==1) + var(cc2$ppvtr.36[cc2$treat==0 & cc2$lths==1]/sum(cc2$treat==0 & cc2$lths==1)))
tes.ns[2,]<-c(mean(cc2$ppvtr.36[cc2$treat==1 & cc2$hs==1]) - mean(cc2$ppvtr.36[cc2$treat==0 & cc2$hs==1]),sum(cc2$treat==1 & cc2$hs==1),sum(cc2$hs==1),var(cc2$ppvtr.36[cc2$treat==1 & cc2$hs==1])/sum(cc2$treat==1 & cc2$hs==1) + var(cc2$ppvtr.36[cc2$treat==0 & cc2$hs==1]/sum(cc2$treat==0 & cc2$hs==1)))
tes.ns[3,]<-c(mean(cc2$ppvtr.36[cc2$treat==1 & cc2$ltcoll==1]) - mean(cc2$ppvtr.36[cc2$treat==0 & cc2$ltcoll==1]),sum(cc2$treat==1 & cc2$ltcoll==1),sum(cc2$ltcoll==1),var(cc2$ppvtr.36[cc2$treat==1 & cc2$ltcoll==1])/sum(cc2$treat==1 & cc2$ltcoll==1) + var(cc2$ppvtr.36[cc2$treat==0 & cc2$ltcoll==1]/sum(cc2$treat==0 & cc2$ltcoll==1)))
tes.ns[4,]<-c(mean(cc2$ppvtr.36[cc2$treat==1 & cc2$college==1]) - mean(cc2$ppvtr.36[cc2$treat==0 & cc2$college==1]),sum(cc2$treat==1 & cc2$college==1),sum(cc2$college==1),var(cc2$ppvtr.36[cc2$treat==1 & cc2$college==1])/sum(cc2$treat==1 & cc2$college==1) + var(cc2$ppvtr.36[cc2$treat==0 & cc2$college==1]/sum(cc2$treat==0 & cc2$college==1)))

round(tes.ns,2)


> round(tes.ns,2)
     [,1] [,2] [,3] [,4]
[1,] 9.30  126 1358 1.80
[2,] 4.06   82 1820 3.31
[3,] 7.87   48  837 5.33
[4,] 4.62   34  366 4.58
>

#these are variances not standard errors in table
#> sqrt(tes.ns[,4])
#[1] 1.341866 1.818015 2.307679 2.140321

## now overall te.s and vars
c(tes.ns[,1])%*%tes.ns[,2]/sum(tes.ns[,2])
        [,1]
[1,] 7.03

c(tes.ns[,1])%*%tes.ns[,3]/sum(tes.ns[,3])
         [,1]
[1,] 6.46

sqrt((c(tes.ns[,2]^2)%*%(tes.ns[,4]))/(sum(tes.ns[,2])^2))
# .90

#num=(126^2)*(1.34^2) + (82^2)*(1.82^2) + (48^2)*(2.31^2) + (34^2)*(2.14^2)
num=(126^2)*(1.80) + (82^2)*(3.31) + (48^2)*(5.33) + (34^2)*(4.58)
denom= (126+82+48+34)^2

sqrt(num/denom)

reg.strat <- lm.all(ppvtr.36 ~ treat + as.factor(educ),data=cc2)
# 6.99


### by birthweight group

tes.ns2=matrix(0,2,3)
tes.ns2[1,]=c(mean(cc2$ppvtr.36[cc2$treat==1 & cc2$bwg==1]) - mean(cc2$ppvtr.36[cc2$treat==0 & cc2$bwg==1]),sum(cc2$treat==1 & cc2$bwg==1),sum(cc2$bwg==1))
tes.ns2[2,]=c(mean(cc2$ppvtr.36[cc2$treat==1 & cc2$bwg==0]) - mean(cc2$ppvtr.36[cc2$treat==0 & cc2$bwg==0]),sum(cc2$treat==1 & cc2$bwg==0),sum(cc2$bwg==0))

         [,1] [,2] [,3]
[1,]  7.63621  142 4159
[2,] 10.39012  235  352

c(tes.ns2[,1])%*%tes.ns2[,2]/sum(tes.ns2[,2])
# 9.28

c(tes.ns2[,1])%*%tes.ns2[,3]/sum(tes.ns2[,3])
# 7.8

####### NOW P-SCORE MATCHING

covs <- c("bw","bwg","hispanic","black","white","b.marr","lths","hs","ltcoll","college","work.dur","prenatal","momage","sex","first","preterm","age","dayskidh" )


######  trying to fit a better model to deal with these issues

### transformed variables now created at top of file

#ps.fit.2 <- glm(treat ~ bwg*as.factor(educ) + as.factor(ethnic)*b.marr +
#                  work.dur + prenatal + preterm + age + momage + sex + 
#                  first + bw + dayskidT + preterm + pretermT + momage + 
#                  momageT + black*(bw + preterm +dayskidT) +  
#                  b.marr*(bw + preterm +dayskidT),
#                data=cc2,family=binomial) 

ps.fit.2 <- glm(treat ~ work.dur + prenatal + preterm + age + momage + sex + 
                  first + bwT + dayskidT + preterm + pretermT + momage +
                  as.factor(ethnic)*b.marr + 
                  b.marr*hs + b.marr*college +
                  b.marr*dayskidT,
                  data=cc2,family=binomial)                 

pscores2 <- ps.fit.2$linear
matches2 <- matching(z=cc2$treat,score=pscores2)
matched2 <- cc2[matches2$matched,]
ps2.mat <- pscores2[matches2$matched]

sum(pscores2[cc2$treat==1] > max(pscores2[cc2$treat==0]))
# 16

reg.ps <- lm(ppvtr.36 ~ treat + bw + bwg + hispanic + black + b.marr + lths + hs + ltcoll + 
               work.dur + prenatal + momage + sex + first + preterm + age + 
               dayskidh, data=matched2)
summary(reg.ps)$coef[2,]
# 10.3 (1.55)

diff.means.matched2=matrix(0,length(covs2),6)
for(i in 1:length(covs2)){
  diff.means.matched2[i,1:2] <- c(mean(matched2[matched2$treat==1,covs2[i]]),mean(matched2[matched2$treat==0,covs2[i]]))
  diff.means.matched2[i,3] <- diff.means.matched2[i,1]-diff.means.matched2[i,2]
  # note in next we still divide it by the standard deviation from full dataset
  diff.means.matched2[i,5] <- sqrt(var(matched2[matched2$treat==1,covs2[i]])/sum(cc2$treat==1) + var(cc2[cc2$treat==0,covs2[i]])/sum(cc2$treat==0))
  diff.means.matched2[i,6] <- sqrt((var(cc2[cc2$treat==1,covs2[i]])+var(cc2[cc2$treat==0,covs2[i]]))/2)
  diff.means.matched2[i,4] <- diff.means.matched2[i,3]/diff.means.matched2[i,6]
}
dimnames(diff.means.matched2) <- list(covs2,c("treat","control","diff","diff.std","se","sd"))
round(diff.means.matched2,2)

covs2C <- c("neg.bw","preterm","dayskidh","age","momage") 

ratio.sds2=matrix(0,length(covs2C),2)
for(i in 1:length(covs2C)){
  ratio.sds2[i,1] <- c(sd(cc2[cc2$treat==1,covs2C[i]])/sd(cc2[cc2$treat==0,covs2C[i]]))
  ratio.sds2[i,2] <- c(sd(matched2[matched2$treat==1,covs2C[i]])/sd(matched2[matched2$treat==0,covs2C[i]]))
}
dimnames(ratio.sds2) <- list(covs2C,c("unmatched ratio","matched ratio"))
round(ratio.sds2,2)


# but we won't use this in the paper -- will use the version in which the pscores
# also had the state indicators
reg.st.ps <- lm(ppvtr.36 ~ treat + hispanic + black + b.marr + lths + hs + ltcoll + 
                  work.dur + prenatal + momage + sex + first + preterm + age + dayskidh + 
                  bw + unemp.rt + st5 + st9 + st12 + st25 + st36 + st42 + st48 + st53,
                  data=matched2)
summary(reg.st.ps)$coef[2,]
# 10.3 (2.26)
########

name=cov.nms
est=diff.means[,3]
sd=diff.means[,6]
est2=diff.means.matched2[,3]
sd2=diff.means.matched2[,6]

J <- 20
#  data.locations = c(1.8,2.8,3.8,4.8,5.8,6.8,9:19)
data.locations = c(1.8,2.8,3.8,4.8,5.8,6.8, 8.9,9.9,10.9,12.2,13.5,14.5,15.5,16.5,17.9,18.9,19.9)
name.range <- .8
x.range <- range (c(est,est2)/c(sd,sd2))
x.range[2] <- x.range[2] +.3
A <- -x.range[1]/(x.range[2]-x.range[1])
B <- 1/(x.range[2]-x.range[1])
height <- .35*J
width <- 3*(name.range+1 )
#

postscript ("balance.both.pretty.ps", horizontal=T, height=4, width=7.5)
par (mar=c(0,0,0,0))
plot (c(-name.range,1), c(3,-J-2), bty="n", xlab="", ylab="",
      xaxt="n", yaxt="n", xaxs="i", yaxs="i", type="n")
text (-name.range, 2, "Predictor", adj=0, cex=1)
text (.5, 2, "Standardized Difference in Means", adj=.5, cex=1)
lines (c(0,1), c(0,0))
lines (c(A,A), c(0,-J-1), lty=2, lwd=.5)
ax <- pretty (x.range)
ax <- ax[(A+B*ax)>0 & (A+B*ax)<1]
segments (A + B*ax, -.1, A + B*ax, .1, lwd=.5)
text (A + B*ax, .7, ax, cex=1)
text (-name.range, -.9, "Child", adj=0, cex=1)
text (-name.range, -8, "Mother", adj=0, cex=1)
text (-name.range+.05, -data.locations, name, adj=0, cex=1)
points (A + B*(est/sd), -data.locations)
points (A + B*(est2/sd2), -data.locations,pch=20, cex=1.5)
#  if (bottom){
#    lines (c(0,1), c(-J-1,-J-1))
#    segments (A + B*ax, -J-1-.1, A + B*ax, -J-1+.1, lwd=.5)
#    text (A + B*ax, -J-1-.7, ax, cex=1)
#  }
dev.off()

##########
### plot of overlap in pscores

##### THIS NEXT IN PAPER (AFTER CONVERTING TO PDF AND CROPPING)
#postscript("ps.overlap.dens.ps", height=3.5, width=5, horizontal=T)
postscript("ihdp.ps.overlap.dens.ps", height=3.5, horizontal=T)
par(mfrow=c(1,2))
par(mar=c(4,2,4,2))
hist(pscores2[cc2$treat==1],xlim=c(-50,max(pscores2)+5),
     main="before matching",xlab="propensity score",freq=F,
     cex.main=1.2,cex.lab=1.2)
hist(pscores2[cc2$treat==0 & pscores2>(-50)],main="",col="grey",border="grey",
     freq=F,add=T)
hist(pscores2[cc2$treat==1],xlim=c(-50,max(pscores2)+5),
     main="before matching",xlab="propensity score",freq=F,
     cex.main=1.2,cex.lab=1.2,add=TRUE)
hist(ps2.mat[matched2$treat==1],xlim=c(min(ps2.mat)-2,7),main="after matching",
     xlab="propensity score",freq=F,cex.main=1.2,cex.lab=1.2)
#hist(ps2.mat[matched2$treat==1],xlim=c(-40,7),main="after matching",
#     xlab="propensity score",freq=F,cex.main=1.2,cex.lab=1.2)
hist(ps2.mat[matched2$treat==0],main="",freq=F,add=T,col="grey",border="grey")
hist(ps2.mat[matched2$treat==1],freq=F,cex.main=1.2,cex.lab=1.2,add=TRUE)
dev.off()

sum(pscores2<(-50))
# 6 controls deleted from the plot to aid in seeing differences across groups
### plot of overlap in pscores FREQUENCY NOT DENSITY

##postscript("ps.overlap.dens.ps", height=3.5, width=5, horizontal=T)
postscript("ihdp.ps.overlap.freq.ps", height=3.5, width=7, horizontal=T)
#pdf("ps.overlap.freq.pdf", height=3.5, width=3.5)
par(mfrow=c(1,2))
par(mar=c(4,2,4,2))
hist(pscores2[cc2$treat==0],xlim=c(min(pscores2),max(pscores2)),main="before matching",xlab="propensity score",
     nclass=20,border="darkgrey",cex.main=1.2,cex.lab=1.2,col="darkgrey")
hist(pscores2[cc2$treat==1],add=T)
hist(ps2.mat[matched2$treat==0],xlim=c(min(ps2.mat)-2,max(pscores2)),main="after matching",xlab="propensity score",
     nclass=10,border="grey",cex.main=1.2,cex.lab=1.2,col="grey")
hist(ps2.mat[matched2$treat==1],nclass=20,add=T)
dev.off()

pscores2I=invlogit(pscores2)
postscript("ihdp.psI.overlap.freq.ps", height=3.5, width=7, horizontal=T)
#pdf("ps.overlap.freq.pdf", height=3.5, width=3.5)
par(mfrow=c(1,2))
par(mar=c(4,2,4,2))
hist(pscores2I[cc2$treat==0&pscores2I>.01],xlim=c(0.01,1),main="before matching",xlab="propensity score",
     nclass=50,border="darkgrey",cex.main=1.2,cex.lab=1.2,col="darkgrey")
hist(pscores2I[cc2$treat==1],add=T,nclass=50)
hist(invlogit(ps2.mat[matched2$treat==0]),xlim=c(0.01,1),main="after matching",xlab="propensity score",
     nclass=50,border="darkgrey",cex.main=1.2,cex.lab=1.2,col="darkgrey")
hist(invlogit(ps2.mat[matched2$treat==1]),nclass=50,add=T)
dev.off()

sum(pscores2I<.01)
# removing the 3652 observations from the plot with pscores<.01

sum(pscores2[cc2$treat==1] > max(pscores2[cc2$treat==0]))

#### UNIVARIATE


### overfit diagnostic
set.seed(1234)
cc2.refit = cbind.data.frame(cc2, treat.rand = rbinom(nrow(cc2),1,mean(cc2$treat)))  
ps.refit.2 <- glm(treat.rand ~ work.dur + prenatal + preterm + age + momage + sex + 
                  first + bwT + dayskidT + preterm + pretermT + momage +
                  as.factor(ethnic)*b.marr + 
                  b.marr*hs + b.marr*college +
                  b.marr*dayskidT,
                data=cc2.refit,family=binomial)                 
#ps.refit.2 <- glm(treat ~ work.dur + prenatal + preterm + age + momage + sex + 
#                  first + bwT + dayskidT + preterm + pretermT + momage +
#                  as.factor(ethnic)*b.marr + b.marr*hs + b.marr*college +
#                  b.marr*dayskidT,
#                 data=cc2.refit,family=binomial) 

pscores2.rf = ps.refit.2$linear
pscores2I.rf=invlogit(pscores2.rf)
#postscript("ihdp.ps.overlap.freq.ps", height=3.5, width=7, horizontal=T)
#pdf("ps.overlap.freq.pdf", height=3.5, width=3.5)
par(mfrow=c(2,1))
par(mar=c(4,2,4,2))
hist(pscores2I.rf[cc2.refit$treat.rand==0],xlim=c(0.0,.25),main="before matching",xlab="propensity score",
     nclass=20,border="darkgrey",cex.main=1.2,cex.lab=1.2,col="darkgrey",freq=FALSE)
hist(pscores2I.rf[cc2.refit$treat.rand==1],add=T,nclass=20,freq=FALSE)
hist(pscores2.rf[cc2.refit$treat.rand==0],xlim=c(min(pscores2.rf),max(pscores2.rf)),main="before matching",xlab="propensity score",
     nclass=20,border="darkgrey",cex.main=1.2,cex.lab=1.2,col="darkgrey",freq=FALSE)
hist(pscores2.rf[cc2.refit$treat.rand==1],add=T,nclass=20,freq=FALSE)
#hist(invlogit(ps2.mat[matched2$treat==0]),xlim=c(0.03,1),main="after matching",xlab="propensity score",
#     nclass=100,border="darkgrey",cex.main=1.2,cex.lab=1.2,col="darkgrey")
#hist(invlogit(ps2.mat[matched2$treat==1]),nclass=100,add=T)

hist(pscores2.rf[cc2.refit$treat.rand==0],xlim=c(-4,-1.5),ylim=c(0,2),main="before matching",
     xlab="propensity score",nclass=20,border="darkgrey",cex.main=1.2,cex.lab=1.2,col="darkgrey",freq=FALSE)
hist(pscores2.rf[cc2.refit$treat.rand==1],add=T,nclass=20,freq=FALSE)

# evaluate wrt the randomization distribution
sum(pscores2.rf[cc2.refit$treat.rand==1] > max(pscores2.rf[cc2.refit$treat.rand==0]))

niters=1000
not.comm=rep(0,niters)
for(i in 1:niters){
  cc2.refit.i = cbind.data.frame(cc2, treat.rand = rbinom(nrow(cc2),1,mean(cc2$treat)))  
  ps.refit.2.i <- glm(treat.rand ~ work.dur + prenatal + preterm + age + momage + sex + 
                    first + bwT + dayskidT + preterm + pretermT + momage +
                    as.factor(ethnic)*b.marr + 
                    b.marr*hs + b.marr*college +
                    b.marr*dayskidT,
                    data=cc2.refit.i,family=binomial)
  pscores2.rf.i = ps.refit.2.i$linear
  not.comm[i] = sum(pscores2.rf.i[cc2.refit.i$treat.rand==1] > max(pscores2.rf.i[cc2.refit.i$treat.rand==0]))
  }


#################
#### and a third version to deal with the geographic location  DOESN'T WORK WELL -- WON'T USE
ps.fit.3 <- glm(treat ~ work.dur + prenatal + preterm + age + momage + sex + 
                  first + bwT + dayskidT + preterm + pretermT + momage +
                  as.factor(ethnic)*b.marr + 
                  b.marr*hs + b.marr*college +
                  b.marr*dayskidT +          
                  st5 + st9 + st12 + st25 + st36 + st42 + st48 + st53,
                data=cc2,family=binomial)                 
#ps.fit.3 <- glm(treat ~ work.dur + prenatal + preterm + age + momage + sex + 
#                 first + bwT + dayskidT + preterm + pretermT + momage +
#                 as.factor(ethnic)*b.marr + b.marr*hs + b.marr*college +
#                 b.marr*dayskidT +               
#                 st5 + st9 + st12 + st25 + st36 + st42 + st48 + st53,
#               data=cc2,family=binomial) 
# 10.22

pscores3 <- ps3.mod$linear
matches3 <- matching(z=cc2$treat,score=pscores3)
matched3 <- cc2[matches3$matched,]
ps3.mat <- pscores3[matches3$matched]

reg.st.ps <- lm(ppvtr.36 ~ treat + hispanic + black + b.marr + lths + hs + ltcoll + 
               work.dur + prenatal + momage + sex + first + preterm + age + dayskidh +
               bw + st5 + st9 + st12 + st25 + st36 + st42 + st48 + st53,
               data=matched3)
summary(reg.st.ps)$coef[2,]
# 8.92 (2.09)

#umemp.rt?
diff.means.matched3=matrix(0,length(covs2),6)
for(i in 1:length(covs2)){
  diff.means.matched3[i,1:2] <- c(mean(matched3[matched3$treat==1,covs2[i]]),mean(matched3[matched3$treat==0,covs2[i]]))
  diff.means.matched3[i,3] <- diff.means.matched3[i,1]-diff.means.matched3[i,2]
  # note in next we still divide it by the standard deviation from full dataset
  diff.means.matched3[i,5] <- sqrt(var(matched3[matched3$treat==1,covs2[i]])/sum(cc2$treat==1) + var(cc2[cc2$treat==0,covs2[i]])/sum(cc2$treat==0))
  diff.means.matched3[i,6] <- sqrt((var(cc2[cc2$treat==1,covs2[i]])+var(cc2[cc2$treat==0,covs2[i]]))/2)
  diff.means.matched3[i,4] <- diff.means.matched3[i,3]/diff.means.matched3[i,6]
}
dimnames(diff.means.matched3) <- list(covs2,c("treat","control","diff","diff.std","se","sd"))
round(diff.means.matched3,2)

ratio.sds3=matrix(0,length(covs2),2)
for(i in 1:length(covs2)){
  ratio.sds3[i,1] <- c(sd(cc2[cc2$treat==1,covs2[i]])/sd(cc2[cc2$treat==0,covs2[i]]))
  ratio.sds3[i,2] <- c(sd(matched3[matched3$treat==1,covs2[i]])/sd(matched3[matched3$treat==0,covs2[i]]))
}
dimnames(ratio.sds3) <- list(covs2,c("unmatched ratio","matched ratio"))
round(ratio.sds3,2)
