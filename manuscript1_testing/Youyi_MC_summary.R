##################################################################################################
# Four predictors 
##################################################################################################


#################################
# Monte Carlo distribution of CVAUC est

library(kyotil)
if(!unix()) setwd("D:/gdrive/auc/cv_auc/R")
proj="est"
ss=c("n25","n40","n80","n500","n2000")
sim.setting=ss%.%"_0_5fold"
reses=lapply(sim.setting, function(x)get.sim.res("res_"%.%proj%.%"/"%.%x, verbose=T))

mypdf(file="MC_distr_cvauc", mfrow=c(2,3))
for (i in 1:length(sim.setting)){
    tmp=reses[[i]]
    plot(density(tmp, adjust=.5), main =ss[i], xlim=range(reses[[1]]), xlab="CV-AUC")
    abline(v=0.5, lty=2)
}
dev.off()


###################################
# comparing CV schemes

library(kyotil)
proj="est"
sim.setting=c("n25","n40","n80")
fit.setting=c("5fold","random5fold","LPO")
settings=c(t(outer(sim.setting, fit.setting, FUN=paste, sep="_"))); if(length(fit.setting)==1) {names(settings)=sim.setting;} else names(settings)=settings
reses=lapply (settings, function(sim.setting) get.sim.res("res_"%.%proj%.%"/"%.%sim.setting, verbose=F))

out=sapply (names(settings), function(s) {
    dat=reses[[s]]
    c(summary(dat), sd=sd(dat))
})
out
colnames(out)=rep(fit.setting,3)
mytex(out, file="tables/summary", col.headers = "\\hline\n  &  \\multicolumn{3}{c}{n25}  &  \\multicolumn{3}{c}{n40}  &  \\multicolumn{3}{c}{n80}  \\\\  \n")



#################################
# bootstrap

library(kyotil)
proj="boot"
sim.setting=c("n25","n40","n80","n500", "n2000")
sim.setting=sim.setting%.%"_0"
fit.setting=c("5fold")
settings=c(t(outer(sim.setting, fit.setting, FUN=paste, sep="_"))); if(length(fit.setting)==1) {names(settings)=sim.setting;} else names(settings)=settings
reses=lapply (settings, function(sim.setting) get.sim.res("res_"%.%proj%.%"/"%.%sim.setting, verbose=F))

out=sapply (names(settings), function(s) {
    dat=reses[[s]]
    c(
          mean(apply(dat, 2, function(x) x["basic.lb"]<0.5 & x["basic.ub"]>0.5 ))
        , mean(apply(dat, 2, function(x) x["perc.lb"]<0.5  & x["perc.ub"]>0.5 ))
        , mean(apply(dat, 2, function(x) x["basic.lb.2"]>0.5))
        , mean(apply(dat, 2, function(x) x["perc.lb.2"]>0.5))
    )
})
out
rownames(out)=c("basic cvg","percentile cvg","basic rej prob","percentile rej prob")
out.boot=out
save(out.boot, file="Rdata/4p_boot.Rdata")

#      n25_0  n40_0  n80_0 n500_0 n2000_0
#[1,] 0.7627 0.7658 0.7715 0.7865  0.7834
#[2,] 0.9570 0.9630 0.9659 0.9672  0.9715
#[3,] 0.0568 0.0521 0.0515 0.0485  0.0461
#[4,] 0.0994 0.0920 0.0844 0.0835  0.0796


#################################
# LeDell CI coverage probabilities

library(kyotil)
proj="LeDell"
truth=c("_0"=0.5, "_0.5"=0.6800024, "_0.75"=0.7585232)
sep="_0"
sim.setting=c("n25","n40","n80","n500", "n2000")
sim.setting=sim.setting%.%sep
fit.setting=c("5fold")
settings=c(t(outer(sim.setting, fit.setting, FUN=paste, sep="_"))); if(length(fit.setting)==1) {names(settings)=sim.setting;} else names(settings)=settings
reses=lapply (settings, function(sim.setting) get.sim.res("res_"%.%proj%.%"/"%.%sim.setting, verbose=F))

out=sapply (names(settings), function(s) {
    dat=reses[[s]]
    c(rowMeans(dat)[1] 
        ,LeDell_cvg=mean(apply(dat, 2, function(x) x["lb"]<truth[sep]  & x["ub"]>truth[sep] ))
        ,LeDell_width=median(dat["ub",]-dat["lb",])
        ,LeDell_pow=mean(dat["lb.2",]>0.5) # lb.2 is 90% CI
    )
})
out
#                 n25_0     n40_0     n80_0     n500_0    n2000_0
#est          0.4998864 0.4993684 0.4993020 0.50040320 0.49991375
#LeDell_cvg   0.7883000 0.8104000 0.8280000 0.83410000 0.83440000
#LeDell_width 0.2327841 0.1889906 0.1539888 0.07160133 0.03578887
#LeDell_pow   0.1365000 0.1237000 0.1149000 0.11640000 0.11510000
out.ledell=out
save(out.ledell, file="Rdata/4p_ledell.Rdata")


#################################
# permtuation test

library(kyotil)
proj="perm"
sim.setting=c("n25","n40","n80","n500","n2000")%.%"_0"
fit.setting=c("5fold")
settings=c(t(outer(sim.setting, fit.setting, FUN=paste, sep="_"))); if(length(fit.setting)==1) {names(settings)=sim.setting;} else names(settings)=settings
reses=lapply (settings, function(sim.setting) get.sim.res("res_"%.%proj%.%"/"%.%sim.setting, verbose=T))

out=sapply (names(settings), function(s) {
    dat=reses[[s]]
    mean(dat["p.value",]<0.05)
})
out
out.perm=out
save(out.perm, file="Rdata/4p_perm.Rdata")


# when using 5-fold in both est and perm distributions
#  n25_0   n40_0   n80_0  n500_0 n2000_0
# 0.0463  0.0470  0.0505  0.0527  0.0514

# when est uses 10x5-fold and perm distribution uses 5-fold
# 0.0378  0.0384  0.0379  0.0411  0.0378


# check estimated CVAUC trained by random forest
library(kyotil)
proj="permrf"
sim.setting=c("n25")%.%c("_0", "_0.1", "_0.2", "_0.3")
fit.setting=c("5fold")
settings=c(t(outer(sim.setting, fit.setting, FUN=paste, sep="_"))); if(length(fit.setting)==1) {names(settings)=sim.setting;} else names(settings)=settings
reses=lapply (settings, function(sim.setting) get.sim.res("res_"%.%proj%.%"/"%.%sim.setting, verbose=T))

sapply(reses, function (x) mean(x))
# close to 0.5 


#################################
# make Table 3, the coverage table, and Table 4, the size table

library(kyotil)
load(file="Rdata/4p_boot.Rdata")
load(file="Rdata/4p_ledell.Rdata")
load(file="Rdata/4p_perm.Rdata")

out.boot
out.ledell
out.perm

tab=rbind(Basic=out.boot["basic cvg",], Percentile=out.boot["percentile cvg",], LeDell=out.ledell["LeDell_cvg",]); tab
colnames(tab)=sub("_0","",colnames(tab))
mytex(tab, file="tables/ci_coverage")


tab=rbind(Percentile=out.boot["percentile rej prob",], Permutation=out.perm, LeDell=out.ledell["LeDell_pow",]); tab
colnames(tab)=sub("_0","",colnames(tab))
mytex(tab, file="tables/type1_error")


#######################################################################
# power, 4P


#################################
# comparing 5-fold and 10x5-fold

library(kyotil)
proj="perm"
sim.setting=c("n25")%.%"_0.5"
fit.setting=c("5fold","10x5fold")
settings=c(t(outer(sim.setting, fit.setting, FUN=paste, sep="_"))); if(length(fit.setting)==1) {names(settings)=sim.setting;} else names(settings)=settings
reses=lapply (settings, function(sim.setting) get.sim.res("res_"%.%proj%.%"/"%.%sim.setting, verbose=T))

out=sapply (names(settings), function(s) {
    dat=reses[[s]]
    mean(dat["p.value",]<0.05)
})
out
#   n25_0.5_5fold n25_0.5_10x5fold
#           0.598            0.635




rm(list=ls())
sep="_0.75" # _0.5 _0.75
library(kyotil)
truth=c("_0.5"=0.6800024, "_0.75"=0.7585232)
sim.setting=c("n25","n40","n80","n500","n2000")%.%sep  
fit.setting=c("5fold")
settings=c(t(outer(sim.setting, fit.setting, FUN=paste, sep="_"))); if(length(fit.setting)==1) {names(settings)=sim.setting;} else names(settings)=settings

proj="LeDell"; reses=lapply (settings, function(sim.setting) get.sim.res("res_"%.%proj%.%"/"%.%sim.setting, verbose=F))
out.1=sapply (names(settings), function(s) {
    dat=reses[[s]]
    c(rowMeans(dat)[1] 
        ,LeDell_cvg=mean(apply(dat, 2, function(x) x["lb"]<truth[sep]  & x["ub"]>truth[sep] ))
        ,LeDell_width=median(dat["ub",]-dat["lb",])
        ,LeDell_pow=mean(dat["lb",]>0.5)
    )
})
out.1

# plot MC distribution
mypdf(mfrow=c(2,3), file="MC_distr"%.%sep)
sapply (names(settings), function(s) {
    dat=reses[[s]]
    plot(density(dat["est",]), main=s)
    abline(v=truth[sep])
})
dev.off()

proj="boot"; reses=lapply (settings, function(sim.setting) get.sim.res("res_"%.%proj%.%"/"%.%sim.setting, verbose=F))
out.2=sapply (names(settings), function(s) {
    dat=reses[[s]]
    c(rowMeans(dat)[1], 
          symm_cvg= mean(apply(dat, 2, function(x) x["symm.lb"]< truth[sep]  & x["symm.ub"]> truth[sep] ))
        , perc_cvg= mean(apply(dat, 2, function(x) x["perc.lb"]< truth[sep]  & x["perc.ub"]> truth[sep] ))
        , basic_cvg=mean(apply(dat, 2, function(x) x["basic.lb"]<truth[sep]  & x["basic.ub"]>truth[sep] ))
        , symm_width= median(dat["symm.ub",]-dat["symm.lb",])
        , perc_width= median(dat["perc.ub",]-dat["perc.lb",])
        , basic_width=median(dat["basic.ub",]-dat["basic.lb",])
        , symm_pow=mean(dat["symm.lb",]>0.5)
        , perc_pow=mean(dat["perc.lb",]>0.5)
        , basic_pow=mean(dat["basic.lb",]>0.5)
        )
})
out.2

# glm
proj="glm"; reses=lapply (settings, function(sim.setting) get.sim.res("res_"%.%proj%.%"/"%.%sim.setting, verbose=F))
out.3=sapply (names(settings), function(s) {
    mean(reses[[s]]<0.05)
})

# permauc
proj="permauc"; reses=lapply (settings, function(sim.setting) get.sim.res("res_"%.%proj%.%"/"%.%sim.setting, verbose=F))
out.4=sapply (names(settings), function(s) {
    dat=reses[[s]]
    mean(dat["p.value",]<0.05)
})

# perm
proj="perm"; reses=lapply (settings, function(sim.setting) get.sim.res("res_"%.%proj%.%"/"%.%sim.setting, verbose=F))
out.5=sapply (names(settings), function(s) {
    dat=reses[[s]]
    mean(dat["p.value",]<0.05)
})

tab=rbind(out.1[1:2,], out.2[2:4,], out.1[3,,drop=F], out.2[5:7,])
tab
mytex(tab, file="tables/ci_coverage"%.%sep)

# don't print LeDell power or basic bootstrap CI power because they undercover
tab=rbind(out.2[8:9,], glm.pow=out.3, perm.auc=out.4, perm.cvauc=out.5)
tab
mytex(tab, file="tables/pow_sim_4p"%.%sep)


#################################
# sim MTCT

rm(list=ls())
library(kyotil)
# 
projs=c("BH","perm_min","primary1","primary2","oracle","perm_lgb","perm_lgb1","perm_dl"); names(projs)=projs
seps=c("_0.0","_0.5")
tab=
sapply(projs, function(proj) {
sapply(seps, function(sep) {
    sim.setting="mtct_"%.%c("n80")%.%sep  
    fit.setting=c("5fold")
    settings=c(t(outer(sim.setting, fit.setting, FUN=paste, sep="_"))); if(length(fit.setting)==1) {names(settings)=sim.setting;} else names(settings)=settings
    reses=lapply (settings, function(sim.setting) get.sim.res("res_"%.%proj%.%"/"%.%sim.setting, verbose=T))
    sapply (names(settings), function(s) {
        if (proj=="BH" | proj=="oracle" | proj=="primary1" | proj=="primary2") {
            mean(reses[[s]]<0.05)
        } else if (startsWith(proj,"perm_")){
            mean(reses[[s]]["p.value",]<0.05)
        }
    })
})
})
rownames(tab)=c("size","power")
tab
mytex(tab, file="tables/pow_sim_mtct")


# examine one set of results in more details
res=get.sim.res("res_perm_lgb1/mtct_n80_0.5_5fold", verbose=T)
summary(res["est",])



###################################################################################################
# sim rv144
###################################################################################################


rm(list=ls())
library(kyotil)
#projs=c("BH","perm_min","primary1","primary2","oracle","perm_lgb","perm_lgb1","perm_dl"); names(projs)=projs
#projs=c("perm_rf", "perm_rf2", "perm_rf3"); names(projs)=projs; seps=c("_1")
projs=c("perm_min", "perm_rf", "perm_rf2"); names(projs)=projs; seps=c("_1")
tab=
sapply(projs, function(proj) {
sapply(seps, function(sep) {
    sim.setting="rv144_"%.%c("n40")%.%sep  
    fit.setting=c("5fold")
    settings=c(t(outer(sim.setting, fit.setting, FUN=paste, sep="_"))); if(length(fit.setting)==1) {names(settings)=sim.setting;} else names(settings)=settings
    reses=lapply (settings, function(sim.setting) get.sim.res("res_"%.%proj%.%"/"%.%sim.setting, verbose=T))
    sapply (names(settings), function(s) {
        if (proj=="BH" | proj=="oracle" | proj=="primary1" | proj=="primary2") {
            mean(reses[[s]]<0.05)
        } else if (startsWith(proj,"perm_")){
            mean(reses[[s]]["p.value",]<0.05)
        }
    })
})
})
tab
rownames(tab)=c("size","power")
tab

mytex(tab, file="tables/pow_sim_rv144_rf")


# examine one set of results in more details
res=get.sim.res("res_perm_rf/rv144_n40_1_5fold", verbose=T)
summary(res["est",])

res.1=get.sim.res("res_perm_rf2/rv144_n40_1_5fold", verbose=T)
summary(res["est",])

res.2=get.sim.res("res_perm_rf3/rv144_n40_1_5fold", verbose=T)
summary(res["est",])
#
## with z: 
#perm_rf._1.rv144_n40_1 perm_min._1.rv144_n40_1
#                  0.840                   0.749
#
## without z:
#perm_rf._1.rv144_n40_1 perm_min._1.rv144_n40_1
#                  0.879                   0.753

# two phase
rm(list=ls())
library(kyotil)
projs=c("perm_min", "perm_rf", "perm_rf2"); names(projs)=projs; 
seps=c("_1")
tab=
sapply(projs, function(proj) {
sapply(seps, function(sep) {
    sim.setting="rv144ph2_"%.%c("n10000")%.%sep  
    fit.setting=c("5fold")
    settings=c(t(outer(sim.setting, fit.setting, FUN=paste, sep="_"))); if(length(fit.setting)==1) {names(settings)=sim.setting;} else names(settings)=settings
    reses=lapply (settings, function(sim.setting) get.sim.res("res_"%.%proj%.%"/"%.%sim.setting, verbose=T))
    sapply (names(settings), function(s) {
        if (proj=="BH" | proj=="oracle" | proj=="primary1" | proj=="primary2") {
            mean(reses[[s]]<0.05)
        } else if (startsWith(proj,"perm_")){
            mean(reses[[s]]["p.value",]<0.05)
        }
    })
})
})
tab
rownames(tab)=c("size","power")
tab



# comparing power between cohort and two phase
rm(list=ls())
library(kyotil)
projs=c("perm_min", "perm_rf", "BH"); names(projs)=projs; 
seps=c("_1")
sims=c("rv144_n40","rv144ph2_n10000"); names(sims)=sims
tab=
sapply(sims, function(sim) {
sapply(projs, function(proj) {
sapply(seps, function(sep) {
    sim.setting=sim%.%sep  
    fit.setting=c("5fold")
    settings=c(t(outer(sim.setting, fit.setting, FUN=paste, sep="_"))); if(length(fit.setting)==1) {names(settings)=sim.setting;} else names(settings)=settings
    reses=lapply (settings, function(sim.setting) get.sim.res("res_"%.%proj%.%"/"%.%sim.setting, verbose=T))
    sapply (names(settings), function(s) {
        if (proj=="BH" | proj=="oracle" | proj=="primary1" | proj=="primary2") {
            mean(reses[[s]]<0.05)
        } else if (startsWith(proj,"perm_")){
            mean(reses[[s]]["p.value",]<0.05)
        }
    })
})
})
})
tab


# size
rm(list=ls())
library(kyotil)
projs=c("perm_min", "perm_rf"); names(projs)=projs; 
seps=c("_0")
sims=c("rv144_n40","rv144ph2_n10000"); names(sims)=sims
tab=
sapply(sims, function(sim) {
sapply(projs, function(proj) {
sapply(seps, function(sep) {
    sim.setting=sim%.%sep  
    fit.setting=c("5fold")
    settings=c(t(outer(sim.setting, fit.setting, FUN=paste, sep="_"))); if(length(fit.setting)==1) {names(settings)=sim.setting;} else names(settings)=settings
    reses=lapply (settings, function(sim.setting) get.sim.res("res_"%.%proj%.%"/"%.%sim.setting, verbose=T))
    sapply (names(settings), function(s) {
        if (proj=="BH" | proj=="oracle" | proj=="primary1" | proj=="primary2") {
            mean(reses[[s]]<0.05)
        } else if (startsWith(proj,"perm_")){
            mean(reses[[s]]["p.value",]<0.05)
        }
    })
})
})
})
tab




rm(list=ls())
library(kyotil)

#res.1=get.sim.res("res_perm_min/rv144_n40_1_5fold")
#res.2=get.sim.res("res_perm_rf/rv144_n40_1_5fold")

#res.3=get.sim.res("res_perm_min/rv144ph2_n10000_1_5fold")
#res.4=get.sim.res("res_perm_rf/rv144ph2_n10000_1_5fold")


tmp=merge(as.data.frame(t(res.1)), as.data.frame(t(res.3)), by="row.names")
head(tmp)
mean(tmp[,3]<.05)
mean(tmp[,5]<.05)

cor(tmp[,3], tmp[,5], method="spearman")
table(tmp[,3]<.05, tmp[,5]<.05)
summary(tmp)


res.1=get.sim.res("res_BH/rv144_n40_1_5fold")
res.3=get.sim.res("res_BH/rv144ph2_n10000_1_5fold")
head(cbind(res.1, res.3))


res.1=get.sim.res("res_test1/rv144_n40_1_5fold")
res.3=get.sim.res("res_test1/rv144ph2_n10000_1_5fold")
mean(apply(res.1, 1, function(x)mean(x<.05)) - apply(res.3, 1, function(x)mean(x<.05)) <0)
summary(apply(res.1, 1, function(x)mean(x<.05)))
summary(apply(res.3, 1, function(x)mean(x<.05)))

summary(apply(res.1,2,min))
summary(apply(res.3,2,min))

table(apply(res.1, 2, which.min))
table(apply(res.3, 2, which.min))

mean(res.1[25,]>res.3[25,])
mean(res.1[1,]>res.3[1,])

mean(apply(res.1,2,min)<.05)
mean(apply(res.3,2,min)<.05)

res.1=get.sim.res("res_BH/rv144_n40_1_5fold")
res.3=get.sim.res("res_BH/rv144ph2_n10000_1_5fold")
mean(res.1<.05)
mean(res.3<.05)


res.1=get.sim.res("res_BH/rv144ph2_n10000_0_5fold")
res.3=get.sim.res("res_test1/rv144ph2_n10000_0_5fold")
mean(res.1<.05)
mean(res.3<.05)


res.1=get.sim.res("res_perm_min/rv144ph2_n10000_1_5fold opt1")
res.2=get.sim.res("res_perm_min/rv144ph2_n10000_1_5fold opt2")
res.3=get.sim.res("res_perm_min/rv144ph2_n10000_1_5fold opt3")
mean(res.1[,2]<.05)
mean(res.2[,2]<.05)
mean(res.3[,2]<.05)


res.1=get.sim.res("res_perm_min/rv144ph2_n10000_1_5fold")
res.2=get.sim.res("res_perm_rf/rv144ph2_n10000_1_5fold")
res.3=get.sim.res("res_perm_rf2/rv144ph2_n10000_1_5fold")

mean(res.1[2,]<.05)
mean(res.2[2,]<.05)
mean(res.3[2,]<.05)


res.1=get.sim.res("res_test2/rv144_n40_1_5fold")
res.3=get.sim.res("res_test2/rv144ph2_n10000_1_5fold")

summary(res.1[25,])
summary(res.3[25,])

summary(res.1[1,])
summary(res.3[1,])
