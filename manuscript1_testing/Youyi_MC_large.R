# simulate datasets like MTCT
rm (list=ls())
library(kyotil); library(MASS); library(boot); library(aucm); library(cvAUC); 
#library(lightgbm); # not able to install by following instructions from https://lightgbm.readthedocs.io/en/latest/R/index.html
library(reticulate)
if(!unix()) setwd("D:/gdrive/MachineLearning/sunwoo/code")
source("sunwoo_shared/helper_cvauc.R")
#source_python("dl.py") #deep learning
#source_python("lgb.py") #lightgbm
    
# process arguments
Args <- commandArgs(trailingOnly=TRUE) 
if (length(Args)==0) {
    # proj: boot for bootstrap, est for est only, LeDell for using cvAUC package
    # sim.setting: mtct_n25_0.5, mtct_n80_0.6, rv144_n40_1, rv144_n40_0
    Args=c(batch.size="10",batch.number="1",sim.setting="rv144_n40_1",fit.setting="5fold",proj="perm_rf2")  
}
myprint(Args)
i=0;
i=i+1; batch.size=as.numeric(Args[i])
i=i+1; batch=as.numeric(Args[i])
seeds=1:batch.size+batch.size*(batch-1); names(seeds)=seeds
myprint(batch, batch.size)
i=i+1; sim.setting=Args[i]
tmp=strsplit(sim.setting,"_")[[1]]
sim.model=tmp[1]
sample.size=tmp[2]
beta=as.numeric(tmp[3]) 
i=i+1; fit.setting=Args[i]
cv.scheme=fit.setting
i=i+1; proj=Args[i]
    
nperm=1e3 # permutation replicates
verbose=ifelse(unix(),0,2)
seed=1 # temp seed for debugging use, does not matter
make.plot=FALSE
if(make.plot) {seeds=1:9; mypdf(mfrow=c(3,3), file=paste0(proj,"_distr_",sim.setting))}
myprint(sample.size, sim.model)

begin=Sys.time()
res=sapply(seeds, simplify="array", function (seed) {
    
    myprint(seed)
    t.0=Sys.time()   
    
    if (sample.size=="n25") {
        n1=25; n0=125
    } else if (sample.size=="n40") {
        n1=40; n0=200
    } else if (sample.size=="n80") {
        n1=80; n0=160
    } else if (sample.size=="n500") {
        n1=500; n0=500
    } else if (sample.size=="n2000") {
        n1=2000; n0=2000
    } else stop("wrong simple.size")
    
    if (sim.model=="mtct") {
        dat=sim.mtct(n1, n0, seed, beta, beta.z.1=1, beta.z.2=1)
    } else if (sim.model=="rv144") {
        if(beta==1) {
            dat=sim.rv144(n1, n0, seed, alpha=-1.9, betas=c(0.9,0.7,1.2,-1,-1.2), beta.z.1=0.5, beta.z.2=0)
        } else if (beta==0) {
            dat=sim.rv144(n1, n0, seed, alpha=-1.9, betas=c(0,0,0,0,0,0), beta.z.1=0, beta.z.2=0)
        }
    } else stop ("wrong sim.model")
    
    # actual number of cases and controls
    n1=sum(dat$y)
    n0=nrow(dat)-n1
    p=ncol(dat)-3 # y, z1, z2
        
    # perform testing
    if (proj=="BH") {
        pvals=sapply (1:p, function(i) last(summary(glm(as.formula("y~z1+z2+x"%.%i), dat, family=binomial))$coef)[4])
        p.adj=p.adjust(pvals, method="BH")
        out=min(p.adj)        
      
    } else if (proj=="oracle") {
#        dat$x.tmp=with(dat, x1+x2+x3+x4+x5)
#        out=summary(glm(y~z1+z2+x.tmp, dat, family=binomial))$coef[2,4]
        fit=glm(y~z1+z2+x1+x2+x3+x4+x5, dat, family=binomial)    # two relevant predictors plus four irrelvant predictors
        fit.0=glm(y~z1+z2+1, dat, family=binomial)    
        out=anova(fit.0, fit, test="Chi")$P[2]
      
    } else if (proj=="primary1") {
        fit=glm(y~z1+z2+x1+x46+x47+x48+x49+x50, dat, family=binomial)    # two relevant predictors plus four irrelvant predictors
        fit.0=glm(y~z1+z2+1, dat, family=binomial)    
        out=anova(fit.0, fit, test="Chi")$P[2]
      
    } else if (proj=="primary2") {
        fit=glm(y~z1+z2+x1+x2+x47+x48+x49+x50, dat, family=binomial)     # two relevant predictors plus four irrelvant predictors
        fit.0=glm(y~z1+z2+1, dat, family=binomial)    
        out=anova(fit.0, fit, test="Chi")$P[2]
      
    } else if (startsWith(proj,"perm")) {   
    # permutation-based tests
    
        # proj-specific implementation
        # dat.b has dat.b$case and dat.b$control
        do.est=function(dat.b){
            X=as.matrix(rbind(dat.b$case, dat.b$control))
            y=c(rep(1,nrow(dat.b$case)), rep(0,nrow(dat.b$control))) # for run_lgb, y needs to be a vector; for run_dl, y needs to be a matrix
            
            if(proj=="perm_min") {
                dat.train=rbind(data.frame(y=1,dat.b$case), data.frame(y=0,dat.b$control))
                pvals=sapply (1:p, function(i) last(summary(glm(as.formula("y~z1+z2+x"%.%i), dat.train, family=binomial))$coef))
                min(pvals)
                
            } else if(startsWith(proj,"perm_rf")) {
            
                splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)                
                
                # no screening, no weighting
                cv.aucs <-  sapply( splits, function(split){
                    dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
                    dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
                    set.seed(123)
                    fit.rf <- randomForest( factor(Y)~., dat.train )
                    pred.rf <- predict( fit.rf, newdata=dat.test, type="prob" )
                    fast.auc( pred.rf[,2], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
                })
                ret=mean(cv.aucs)
                
                if (proj=="perm_rf2") {
                    # screening, no weighting
                    cv.aucs.2 <-  sapply( splits, function(split){
                        dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
                        dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
                        set.seed(123)
                        # screening
                        predictors=setdiff(names(dat.train),"Y")
                        pvals=sapply(predictors, function (x) last(c(summary(glm(as.formula("Y~"%.%x), dat.train, family=binomial()))$coef)) )
                        predictors=predictors[pvals < .1]
                        fit.rf <- randomForest( as.formula("factor(Y)~"%.% concatList(predictors,"+")), dat.train )
                        pred.rf <- predict( fit.rf, newdata=dat.test, type="prob" )
                        fast.auc( pred.rf[,2], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
                    })
                    ret.2=mean(cv.aucs.2)
                    # maximum of with and without screening
                    ret=max(ret, ret.2)              
                }
                ret

            } else if(startsWith(proj,"perm_dl")) {
                # deep learning through the python script
                if (sim.model=="mtct") {
                    if (proj=="perm_dl") {
                        k1=10; k2=3
                    }
                } else if (sim.model=="rv144") {
                    if (proj=="perm_dl") {
                        #k1=10; k2=3# power 0.05
                        #k1=100; k2=10# power 0.05
                        k1=20; k2=5# power 0.05
                    }
                }
                out=run_dl(X, as.matrix(y), as.integer(k1), as.integer(k2), seed=0, verbose=verbose)
                max(unlist(out)) # cv-auc
                
            } else if(startsWith(proj,"perm_lgb")) {
                # there are two ways to run lgb
                # (1):  call a python script
                depth=3; leaves=5
#                res <- run_lgb(X, y, depth, leaves, seed)
                
                # (2) through lightgbm R library
                if (proj=="perm_lgb") {
                    lgb.grid = list(objective = "binary",
                                    num_leaves=5, # max number of leaves
                                    min_data_in_leaf=10, 
                                    max_depth=3, 
                                    learning_rate=0.05,
                                    is_unbalance = TRUE,
                                    metric = "auc")
                } else if (proj=="perm_lgb1") {
                    lgb.grid = list(objective = "binary",
                                    num_leaves=5, # max number of leaves
                                    min_data_in_leaf=10, 
                                    max_depth=5, 
                                    learning_rate=0.05,
                                    is_unbalance = TRUE,
                                    metric = "auc")
                } else stop("wrong proj: "%.%proj)
                
                dtrain <- lgb.Dataset(X, label=y)
                
                save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
                if (class(save.seed)=="try-error") { set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }                        
                set.seed(1)# for lgb.cv
                
                # verbose -1 suppress the warnings
                # num_threads 1 makes it much faster on Linux for some reason!
                # nrounds is number of trees
                lgb.model.cv <- lgb.cv(lgb.grid, dtrain, verbose=-1, num_threads=1
                    , nrounds=1000, early_stopping_rounds=5, nfold=5, stratified=T) 
                best.iter = lgb.model.cv$best_iter
                print(best.iter)
                lgb.model.cv$best_score
                
                assign(".Random.seed", save.seed, .GlobalEnv)     
                                
                # need to define dtrain again if run on windows, otherwise there is error
                # redefine this make it run twice as fast on Linux!
                dtrain <- lgb.Dataset(X, label=y)
                lgb.model <- lgb.train(lgb.grid, dtrain, verbose=-1, num_threads=1
                    , nrounds=best.iter)
                fast.auc(predict(lgb.model,as.matrix(rbind(dat.b$case, dat.b$control))), c(rep(1,n1), rep(0,n0)), reverse.sign.if.nece = F)
    
            } else stop("wrong proj")
        }
        
        dat.tmp=list(case=subset(dat, y==1, select=-y), control=subset(dat, y==0, select=-y))
        dat.case.control=rbind(dat.tmp$case, dat.tmp$control)
        est=do.est(dat.tmp)        
        
        N=n0+n1        
        if (choose(N,n1)<=nperm) p.method="exact" else p.method="Monte Carlo"
        if(p.method=="exact"){        
            indexes=combn(1:N, n1) 
            ref.distr=apply(indexes, 2, function(ind.1) {
                ind.0=setdiff(1:N, ind.1)
                dat.b=list()
                dat.b$case=   cbind(dat.case.control[1:n1,c("z1","z2")],    subset(dat.case.control[ind.1,], select=-c(z1,z2)) )
                dat.b$control=cbind(dat.case.control[1:n0+n1,c("z1","z2")], subset(dat.case.control[ind.0,], select=-c(z1,z2)) )
#                dat.b$case=   dat.case.control[ind.1,]
#                dat.b$control=dat.case.control[ind.0,]
                do.est(dat.b)
            })
            
        } else if (p.method=="Monte Carlo") {
            #### save rng state before set.seed in order to restore before exiting this function
            save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
            if (class(save.seed)=="try-error") { set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }                        
            set.seed(1)
            
            ref.distr=replicate(nperm, {
                ind.1=sample(1:N, n1)    # cases in the permutated dataset
                ind.0=setdiff(1:N,ind.1) # controls in the permutated dataset
                dat.b=list()
                dat.b$case=   cbind(dat.case.control[1:n1,c("z1","z2")],    subset(dat.case.control[ind.1,], select=-c(z1,z2)) )
                dat.b$control=cbind(dat.case.control[1:n0+n1,c("z1","z2")], subset(dat.case.control[ind.0,], select=-c(z1,z2)) )
                do.est(dat.b)
            })            
                
            #### restore rng state 
            assign(".Random.seed", save.seed, .GlobalEnv)     
        } else stop("wrong p.method")
        
        # plot permutation samples distribution
        if(make.plot) {plot(density(ref.distr), main=seed); abline(v=est, lty=2); abline(v=0.5)}
        
        # get pvalue
        if(proj=="perm_min") {
            p.value=mean(ref.distr<est)
        } else {
            #p.value=ifelse(est>0.5, mean(ref.distr>est)+mean(ref.distr<1-est), mean(ref.distr<est)+mean(ref.distr>1-est))            
            p.value = mean(ref.distr>=est)
        } 
        
        out=c(est=est, p.value=p.value)
    
    } else stop("wrong proj")
    
    gc()
    if(verbose) print("time used: "%.%format(Sys.time()-t.0))      
        
    out
})
if(make.plot) dev.off()
print(date()%.%". Time used: "%.%format(Sys.time()-begin))

# save res
foldername="res_"%.%proj%.%"/"; if(!file.exists(foldername)) dir.create(foldername)
foldername=foldername%.%sim.setting%.%"_"%.%fit.setting%.%"/"; if(!file.exists(foldername)) dir.create(foldername)
save (res, file=foldername%.%"/batch"%.%formatInt(batch, 3)%.%".Rdata")
