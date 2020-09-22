rm (list=ls()) 
library(kyotil) 
library(MASS) 
library(boot) 
library(aucm)
library(cvAUC)
library(MLmetrics)
#library(lightgbm)
library(reticulate)
if (unix()) {
    source("../helper_cvauc.R")
    source("../helper_rf_hvtn.R")
} else {
    setwd("D:/gdrive/MachineLearning/sunwoo/code")
    source("sunwoo_shared/helper_cvauc.R")
    source("sunwoo_shared/helper_rf_hvtn.R")
}



# make a plot of correlation between predictors for an RV144 dataset
dat=sim.rv144(40, 200, seed=1, alpha=-1.9, betas=c(0.9,0.7,1.2,-1,-1.2), beta.z.1=0.5, beta.z.2=0)        
cor=cor(dat[,-(1:3)])

# this requires a hack in DMHeatMap to create the vertical lines
myfigure()
    hU=DMHeatMap(cor, trace="none", symm=T,dendrogram="none", col=RColorBrewer::brewer.pal(length(breaks)-1,"RdYlGn"), distfun = function(c) as.dist(1 - c), axis=F,  breaks=c(-1,-.7,-.5,-.3,-.1,.1,.3,.5,.7,1), margins=c(0,0), key = F, Rowv=NA, lower.left.only=FALSE, heatmapOnly=T) 
mydev.off(file="../writing/sunwoo_manus1/input/rv144_cor")


# the following does not work because abline works not intuitively
heatmap(cor, symm=T, col=RColorBrewer::brewer.pal(length(breaks)-1,"RdYlGn"), distfun = function(c) as.dist(1 - c), breaks=c(-1,-.7,-.5,-.3,-.1,.1,.3,.5,.7,1), margins=c(0,0), Rowv=NA) 
abline(v=(cumsum(c(7,2,15,10))+.5)/150)



##
beta.x=1 # 1 , .3
beta.zs=c(1,.5,0)
out=sapply(beta.zs, function(beta.z) {
    alpha=1
    n=200
    res=(sapply(1:100, function(seed){        
        set.seed(seed)
        z=rnorm(n)
        x=rnorm(n)
        eta=alpha+z*beta.z+x*beta.x
        y=rbern(n, expit(eta))
        mean(y)
        dat=data.frame(z,x,y)
        pval=summary(glm(y~z+x,dat,family=binomial()))$coef[3,4]
        c(z=fast.auc(z, y), x=fast.auc(x, y), eta=fast.auc(eta, y), p=pval)
    }))
    apply(res, 1, median)
})
colnames(out)="beta.z="%.%beta.zs
out





# get true AUC

# sep=0.5
n=1e6
seed=3
dat=sim.1(n, n, seed=seed, sep=.5)
dat.train=rbind(data.frame(Y=1,dat$case),   data.frame(Y=0,dat$control))
dat=sim.1(n, n, seed=987654+seed, sep=.5)
dat.test =rbind(data.frame(Y=1,dat$case),   data.frame(Y=0,dat$control))
fit=glm(Y~X1+X2+X3+X4, dat.train, family=binomial)    
fast.auc(predict(fit, newdata=dat.test), dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE) 
# seed 1: 0.6799800
# seed 2: 0.6799464
# seed 3: 0.6800809
# average: 0.6800024


# sep=0.75
n=1e6
seed=3
sep=0.75
dat=sim.1(n, n, seed=seed, sep=sep)
dat.train=rbind(data.frame(Y=1,dat$case),   data.frame(Y=0,dat$control))
dat=sim.1(n, n, seed=987654+seed, sep=sep)
dat.test =rbind(data.frame(Y=1,dat$case),   data.frame(Y=0,dat$control))
fit=glm(Y~X1+X2+X3+X4, dat.train, family=binomial)    
fast.auc(predict(fit, newdata=dat.test), dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE) 
# seed 1: 0.7585258
# seed 2: 0.7584554
# seed 3: 0.7585883
# average: 0.7585232


# plot sep
dat=sim.1(n1=80, n0=160, seed=1, sep=0.5)
dat.train=rbind(data.frame(Y=1,dat$case),   data.frame(Y=0,dat$control))
with(dat.train, plot(X2~X1, col=Y+1))
myboxplot(X2~Y, dat.train)
fit=glm(Y~X1+X2+X3+X4, dat.train, family=binomial)    
summary(fit)
fit.0=glm(Y~1, dat.train, family=binomial)    
anova(fit.0, fit, test="Chi")$P[2]


###############################################################################################
# MTCT dataset

library(MTCT)
# list of all immune variables
assays=c("V3_", "NAb", "IgG", "Avdt_", "ADCC", "IgA"); names(assays)=assays # CD4 not included to mimic the mtct pre-specified plan, but CD4 also highly correlated with V3 and NAb_tier 1
all.vars=unlist(sapply (assays, function (label) {
    vv=get.var(full.m, label)
    vv=vv[!endsWith(vv, "_cat") & !endsWith(vv, "_score") & !startsWith(vv, "V3_BioV3B_") & !startsWith(vv, "V3_BioV3Bnew_") & !startsWith(vv, "V3_BioV3M_") & !startsWith(vv, "Avdt_lru") & !startsWith(vv, "Avdt_lkoff") & !contain(vv, ".imp")]
}))
n.vars=length(all.vars)
# remove some variables including controls and new V3B etc
all.vars.2=setdiff(all.vars, c("NAb_SVAMLV","IgG_MuLVgp70_His6","V3_BioV3Bnew","V3_gp70MNV3_auc","V3_BioV3Bnew_ce"))
all.vars.2

library(immuCorr)
cor.heatmap(full.m[,all.vars.2])



dat=sim.mtct(80,160,1,0.7)

fast.auc(dat$x1, dat$y) #0.587

fit=glm(y~x1, dat, family=binomial); summary(fit)$coef
fit=glm(y~x2, dat, family=binomial); summary(fit)$coef
fit=glm(y~x3, dat, family=binomial); summary(fit)$coef
fit=glm(y~x4, dat, family=binomial); summary(fit)$coef
fit=glm(y~x5, dat, family=binomial); summary(fit)$coef
fit=glm(y~x6, dat, family=binomial); summary(fit)$coef


# other parameters
#                min_sum_hessian_in_leaf = 1,
#                feature_fraction = 0.7,
#                bagging_fraction = 0.7,
#                bagging_freq = 5,
#                min_data = 100,
#                max_bin = 50,
#                lambda_l1 = 8,
#                lambda_l2 = 1.3,
#                min_data_in_bin=100,
#                min_gain_to_split = 10,

lgb.grid.1 = list(objective = "binary",
                max_depth=3, 
                num_leaves=5, # max number of leaves
                min_data_in_leaf=10, 
                learning_rate=0.05,
                is_unbalance = TRUE,
                metric = "auc"
)

lgb.grid.2 = list(objective = "binary",
                max_depth=5, # 5 is best judging by validation 
                num_leaves=10, # max number of leaves
                min_data_in_leaf=10, 
                learning_rate=0.05,
                is_unbalance = TRUE,
                metric = "auc"
)

lgb.grid=lgb.grid.2
rand.y=(dat$y)
set.seed(1) # for cv
dtrain <- lgb.Dataset(as.matrix(dat[,-1]), label=rand.y)
lgb.model.cv <- lgb.cv(lgb.grid, dtrain, nrounds=1000, early_stopping_rounds=5, nfold=5, stratified=T) # nrounds is number of trees
best.iter = lgb.model.cv$best_iter; print(best.iter)

#
dtrain <- lgb.Dataset(as.matrix(dat[,-1]), label=rand.y)
lgb.model <- lgb.train(lgb.grid, dtrain, nrounds=best.iter)
#
fast.auc(predict(lgb.model,as.matrix(dat[,-1])), rand.y, reverse.sign.if.nece = F)
lgb.importance(lgb.model, percentage = F)


# can we use different learner?


library(rjson)
library(data.tree)
json_model =lgb.dump(lgb.model)
result <- fromJSON(json_model)
#write(json_model, file="json_model.txt")
node=as.Node(result)
print(node, "split_feature")


##############
# integrate python

library(reticulate)
py_available(TRUE)

os <- import("os")
setwd("D:/gdrive/auc/cv_auc/R")
os$getcwd()

dat=sim.mtct(80,160,1,0.7)

mywrite.csv(dat[,1:4],file="data/example")

fit=glm(y~x1+x2+x3, dat, family=binomial)

source_python("lgb.py")

res=
sapply(2:5, function(depth) {
sapply(2:5, function(leaves) {
    mean(sapply(1:10, function(seed) last(run_lgb(subset(dat, select=-y), subset(dat, select=y), depth, leaves, seed=seed)))  )
})
})
res
#          [,1]      [,2]      [,3]      [,4]
#[1,] 0.5518124 0.5518124 0.5518124 0.5518124
#[2,] 0.5692250 0.5692250 0.5692250 0.5692250
#[3,] 0.5603287 0.5679445 0.5679445 0.5679445
#[4,] 0.5603287 0.5565811 0.5504495 0.5504495
colMeans(res)

depth=3; leaves=3
# do lgb.cv through R. the random number generators are different between R and python, so we don't get the same results as python call
lgb.grid = list(objective = "binary",
                max_depth=depth, # 5 is best judging by validation 
                num_leaves=leaves, # max number of leaves
                min_data_in_leaf=10, 
                learning_rate=0.05,
                is_unbalance = TRUE,
                metric = "auc"
)
rand.y=(dat$y)
set.seed(0) # for cv
dtrain <- lgb.Dataset(as.matrix(dat[,-1]), label=rand.y)
lgb.model.cv <- lgb.cv(lgb.grid, dtrain, nrounds=1000, early_stopping_rounds=10, nfold=5, stratified=T) # nrounds is number of trees
best.iter = lgb.model.cv$best_iter; print(best.iter)


# call dl
library(kyotil) 
library(reticulate)
if(!unix()) setwd("D:/gdrive/auc/cv_auc/R")
source("helper_cvauc.R")
source_python("dl.py")
dat=sim.mtct(80,160,1,0.7)
X=subset(dat, select=-y)
y=subset(dat, select=y)
k1=10; k2=3; seed=0
out=run_dl(as.matrix(X), as.matrix(y), as.integer(k1), as.integer(k2), seed); unlist(out)  #as.integer is needed for type in python




##########################################################################################################
# test sim.rv144
seeds=1:100
betas=c(0.9,0.7,0.7,0.5,-1)

# distributions
dat=sim.rv144(40,200,seed=2, betas)
apply(dat[,-1],2,sd)
round(cor(dat[,-1]),2)[35:64,35:64]


# grp1
res=sapply (seeds, function(seed) {
    dat=sim.rv144(40,200,seed, betas)
    sapply (1:7, function(i) round(summary(glm(as.formula("y~x"%.%i), dat, family=binomial))$coef[2,4],2))
})
rowMeans(res<0.05)

# grp2
res=sapply (seeds, function(seed) {
    dat=sim.rv144(40,200,seed, betas)
    sapply (8:9, function(i) round(summary(glm(as.formula("y~x"%.%i), dat, family=binomial))$coef[2,4],2))
})
rowMeans(res<0.05)

# grp3
res=sapply (seeds, function(seed) {
    dat=sim.rv144(40,200,seed, betas)
    sapply (10:24, function(i) round(summary(glm(as.formula("y~x"%.%i), dat, family=binomial))$coef[2,4],2))
})
rowMeans(res<0.05)
# grp 3 log transformed
res=sapply (seeds, function(seed) {
    dat=sim.rv144(40,200,seed, betas)
    sapply (10:24, function(i) round(summary(glm(as.formula("y~log(x"%.%i%.%")"), dat, family=binomial))$coef[2,4],2))
})
rowMeans(res<0.05)

# grp4
res=sapply (seeds, function(seed) {
    dat=sim.rv144(40,200,seed, betas)
    sapply (25:34, function(i) round(summary(glm(as.formula("y~x"%.%i), dat, family=binomial))$coef[2,4],2))
})
rowMeans(res<0.05)

# grp 2 and 4 interaction
res=sapply (seeds, function(seed) {
    dat=sim.rv144(40,200,seed, betas)
    sapply (25:34, function(i) {
        round(summary(glm(as.formula("y~x8*x"%.%i), dat, family=binomial))$coef[4,4],2)
    })
})
rowMeans(res<0.05)




#################################################
# rv144cc dataset

library(RV144cc)

prim.var=c("ADCC","Avidity","IgA","ICS","V2","NAB")

tmp=names(allcc)
tmp=tmp[!endsWith(tmp, "_cat") & !endsWith(tmp, "_imp")]# toggle_luminex_cytokines_mcelrath_wk26_imp  is binary for some reason

tmp[startsWith(tmp, "secondary")]

markers=c(prim.var, tmp[startsWith(tmp, "toggle")], )

gender + behavioral_risk
# behavioral_risk is a factor with three levels

library(splines)
fit=glm(y~gender + behavioral_risk+ secondary_jm_lum_il10_wk26, allcc, family="binomial"); summary(fit)

dat=subset(allcc, !is.na(secondary_jm_lum_il10_wk26))


fit.1=glm(y~gender + behavioral_risk+IgA+V2, dat, family="binomial"); summary(fit.1)
fit.2=glm(y~gender + behavioral_risk+IgA+ns(V2, df=2), dat, family="binomial"); summary(fit.2)
fit.3=glm(y~gender + behavioral_risk+IgA+ns(V2, df=3), dat, family="binomial"); summary(fit.3)

xx=quantile(dat$V2,seq(.1,.9,length=50))
plot(xx,  predict(fit.3, newdata=data.frame(gender="F", behavioral_risk="Low",IgA=0,V2=xx) ), type="l", ylab="logit(risk)", xlab="V2")
lines(xx, predict(fit.2, newdata=data.frame(gender="F", behavioral_risk="Low",IgA=0,V2=xx) ))

myfigure(mfrow=c(1,2))
    a="V2"; binaryloess(allcc[[a]], allcc$y, span=.7, xlab=a, ylab="infection probability", scale="linear")
    binaryloess(allcc[["secondary_jm_lum_il10_wk26"]], allcc$y, s=.7, xlab="IL10", ylab="infection probability", scale="linear")


myboxplot(secondary_jm_lum_il10_wk26~y, dat)

allcc$V2sq=allcc$V2**2
allcc$V2cu=allcc$V2**3
fit1=run.cc(~V2+V2sq+V2cu, allcc, ret.risk=T)
plot.risk("V2", attr(fit1,"risk"))


allcc$IL10=allcc$secondary_jm_lum_il10_wk26
allcc$IL10sq=allcc$IL10**2
allcc$IL10cu=allcc$IL10**3
fit1=run.cc(~IL10+IL10sq+IL10cu, allcc, ret.risk=T)
plot.risk("IL10", attr(fit1,"risk"))


library(chngpt)
fit=chngptm(y~gender + behavioral_risk, ~secondary_jm_lum_il10_wk26, allcc, type="M10", family="binomial")
fit=chngptm(y~gender + behavioral_risk, ~secondary_jm_lum_il10_wk26, allcc, type="M20", family="binomial")


dat$cut=cut(dat$V2, quantile(dat$V2, seq(0,1,length=5), na.rm=T))
table.prop(dat$y, dat$cut)

a="V2"
dat$cut=cut(dat[[a]], quantile(dat[[a]], seq(0,1,length=6), na.rm=T))
table.prop(dat$y, dat$cut)

fit=glm(y~gender + behavioral_risk+ V2, allcc, family="binomial"); summary(fit)
fit=chngptm(y~gender + behavioral_risk, ~V2, allcc, type="M12c", family="binomial"); summary(fit)


a="IgGFI26bs_gp70V1V2CaseA2"; hist(BAMA.m[[a]])
allcc[[a]]=BAMA.m[[a]][match(allcc$ptid,BAMA.m$ptid)]
fit=glm(y~gender + behavioral_risk+ IgGFI26bs_gp70V1V2CaseA2, allcc, family="binomial"); summary(fit)
fit=chngptm(y~gender + behavioral_risk, ~IgGFI26bs_gp70V1V2CaseA2, allcc, type="M01", family="binomial"); summary(fit)
