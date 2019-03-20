### "truth" for this group is 7.4

cc <- read.table("data/ihdp.nlsy.imp1.txt",header=T,sep="\t")
cc$treat <- cc$treat-1

cc$ethnic = cc$hispanic
cc$ethnic[cc$black==1]=2
cc$ethnic[cc$white==1]=3

cc$educ = cc$lths
cc$educ[cc$hs==1]=2
cc$educ[cc$ltcoll==1]=3
cc$educ[cc$college==1]=4
cc$educ3=cc$educ
cc$educ3[cc$educ>2]=cc$educ3[cc$educ>2]-1

cc$bwg=(cc$bw>2000)*1

cc$state=cc$st5
cc$state[cc$st9==1]=2
cc$state[cc$st12==1]=3
cc$state[cc$st25==1]=4
cc$state[cc$st36==1]=5
cc$state[cc$st42==1]=6
cc$state[cc$st48==1]=7
cc$state[cc$st53==1]=8

cc$state2=cc$state
cc$state2[cc$st5==1]=0

cc$state3=cc$state
cc$state3[cc$st53==1]=0
#covs <- c("hispanic","black","white","b.marr","lths","hs","ltcoll","college","work.dur","prenatal","momage","sex","first","preterm","age","dayskidh","bw","unemp.rt","st5","st9","st12","st25","st36","st42","st48","st53")

# code for looking at nlsy/ihdp high quality child care data
cc2<-cc[cc$bw>1500,]

cc2$neg.bw = 2500 - cc2$bw
cc2$no.prenatal = 1-cc2$prenatal
cc2$b.unmarr = 1-cc2$b.marr

save(cc2,file="data/cc2.Rdata")
