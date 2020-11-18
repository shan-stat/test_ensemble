# simulate datasets like MTCT
rm (list=ls())
library(kyotil); library(MASS); library(boot); library(aucm); library(cvAUC); library(survey); library(ranger)
#library(lightgbm); # not able to install by following instructions from https://lightgbm.readthedocs.io/en/latest/R/index.html
library(reticulate)
if (unix()) {
    source("../helper_cvauc.R")
    source("../helper_rf_hvtn.R")
} else {
    setwd("D:/gdrive/MachineLearning/sunwoo/code")
    source("sunwoo_shared/helper_cvauc.R")
    source("sunwoo_shared/helper_rf_hvtn.R")
}
#source_python("dl.py") #deep learning
#source_python("lgb.py") #lightgbm
    
# process arguments
Args <- commandArgs(trailingOnly=TRUE) 
if (length(Args)==0) {
    # proj: boot for bootstrap, est for est only, LeDell for using cvAUC package
    # sim.setting: mtct_n25_0.5, mtct_n80_0.6, rv144_n40_0, rv144_n40_1, rv144ph2_n10000_1
    Args=c(batch.size="10",batch.number="1",sim.setting="rv144ph2_n10000_1",fit.setting="5fold",proj="test2")  
}
myprint(Args)
i=0;
i=i+1; batch.size=as.numeric(Args[i])
i=i+1; batch=as.numeric(Args[i])
seeds=1:batch.size+batch.size*(batch-1); names(seeds)=seeds
myprint(batch, batch.size)
i=i+1; sim.setting=Args[i]; tmp=strsplit(sim.setting,"_")[[1]]
    sim.model=tmp[1]
    sample.size=tmp[2]
    beta=as.numeric(tmp[3]) 
i=i+1; fit.setting=Args[i]
    cv.scheme=fit.setting
    fit2ph=endsWith(sim.model,"ph2")
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
    } else if (sample.size=="n10000") {
        n1=NULL; n0=NULL; n=1e4
    } else stop("wrong simple.size")
    
    if (sim.model=="mtct") {
        dat=sim.mtct(n1, n0, seed, beta, beta.z.1=1, beta.z.2=1)
    } else if (startsWith(sim.model,"rv144"))  {
        if (sim.model=="rv144") {
            dat=sim.rv144(n1, n0, seed, alpha=-1.9,                     betas=if(beta==1) c(0.9,0.7,1.2,-1,-1.2) else c(0,0,0,0,0,0), beta.z.1=0.5, beta.z.2=0)        
        } else if (sim.model=="rv144ph2") {
            # dat=sim.rv144(n1, n0, seed, alpha=if(beta==1) -7.6 else -6, betas=if(beta==1) c(0.9,0.7,1.2,-1,-1.2) else c(0,0,0,0,0,0), beta.z.1=0.5, beta.z.2=0, n=n)        
            # min rf rf2: 1, 1, 1  
            # dat=sim.rv144(n1, n0, seed, alpha=if(beta==1) -7.6 else -6, betas=if(beta==1) c(0.5,0.5, 0.2,-.5,-1.5)  else c(0,0,0,0,0,0), beta.z.1=0.5, beta.z.2=0, n=n) 
            # min rf rf2: 0.7, 0.7, 0.7  
            # dat=sim.rv144(n1, n0, seed, alpha=if(beta==1) -7.6 else -6, betas=if(beta==1) c(0.5,0.5, 0.15,-.15,-.15)  else c(0,0,0,0,0,0), beta.z.1=0.5, beta.z.2=0, n=n)   
            # min rf rf2: 0.06, 0.03, 0.02  
            # dat=sim.rv144(n1, n0, seed, alpha=if(beta==1) -7.6 else -6, betas=if(beta==1) c(0.5,0.5, 0.5,-.5,-.5)  else c(0,0,0,0,0,0), beta.z.1=0.5, beta.z.2=0, n=n)   
            # min rf rf2: 0.09, 0.09, 0.08  
            # dat=sim.rv144(n1, n0, seed, alpha=if(beta==1) -7.6 else -6, betas=if(beta==1) c(0.5,0.5, 0.5,-.5,-1)  else c(0,0,0,0,0,0), beta.z.1=0.5, beta.z.2=0, n=n)   
            # min rf rf2: 0.19, 0.26, 0.26  
            # dat=sim.rv144(n1, n0, seed, alpha=if(beta==1) -7.6 else -6, betas=if(beta==1) c(0.5,0.5, 0.8,-.8,-1.2)  else c(0,0,0,0,0,0), beta.z.1=0.5, beta.z.2=0, n=n)   
            # min rf rf2: 0.82, 0.76, 0.74  
            # dat=sim.rv144(n1, n0, seed, alpha=if(beta==1) -7.6 else -6, betas=if(beta==1) c(0.5,0.5, 0.6,-.8,-1.5)  else c(0,0,0,0,0,0), beta.z.1=0.5, beta.z.2=0, n=n)   
            # min rf rf2: .98, .97, .99  
            # dat=sim.rv144(n1, n0, seed, alpha=if(beta==1) -7.6 else -6, betas=if(beta==1) c(0.5,0.5, 0.6,-.6,-1.2)  else c(0,0,0,0,0,0), beta.z.1=0.5, beta.z.2=0, n=n) 
            # min rf rf2: .47, .44, .48  
            # dat=sim.rv144(n1, n0, seed, alpha=if(beta==1) -7.6 else -6, betas=if(beta==1) c(0.5,0.5, 0.6,-.6,-1.4)  else c(0,0,0,0,0,0), beta.z.1=0.5, beta.z.2=0, n=n) 
            # min rf rf2: .75, .73, .76  
            #dat=sim.rv144(n1, n0, seed, alpha=if(beta==1) -7.6 else -6, betas=if(beta==1) c(0.5,0.5, 0.6,-.8,-1.4)  else c(0,0,0,0,0,0), beta.z.1=0.5, beta.z.2=0, n=n) 
            # min rf rf2: .93, .90, .93  
            #dat=sim.rv144(n1, n0, seed, alpha=if(beta==1) -7.6 else -6, betas=if(beta==1) c(0.5,0.5, 0.8,-.6,-1.4)  else c(0,0,0,0,0,0), beta.z.1=0.5, beta.z.2=0, n=n) 
            # min rf rf2: .77, .74, .76  
            # dat=sim.rv144(n1, n0, seed, alpha=if(beta==1) -7.6 else -6, betas=if(beta==1) c(0.5,0.5, 0.8,0,-1.5)  else c(0,0,0,0,0,0), beta.z.1=0.5, beta.z.2=0, n=n) 
            # min rf rf2: .18, .32, xx
            dat=sim.rv144(n1, n0, seed, alpha=if(beta==1) -7.6 else -6, betas=if(beta==1) c(0.5,0.5, 0.8,0,-2)  else c(0,0,0,0,0,0), beta.z.1=0.5, beta.z.2=0, n=n) 
            # min rf rf2: .46, .89, .88 
        }
        f = "as.factor(Y)~z1+z2+" %.% concatList(names(dat)[startsWith(names(dat),"x")],"+")    
    } else stop ("wrong sim.model")
    
    # actual number of cases and controls
    N=nrow(dat)
    n1=sum(dat$y)
    n0=N-n1
    p=sum(startsWith(names(dat),"x"))
    myprint(n1,n0,p)


    # perform testing
    if (proj=="BH") {
        pvals=sapply (1:p, function(i) {
            set.seed(123)
            if(!fit2ph) {
                fit=glm(as.formula("y~z1+z2+x"%.%i), dat, family=binomial)
            } else {
                # option 1: use survey package
                #dstrat<-svydesign(id=~1,strata=~bstrat, weights=~wt, data=dat) # stratified sampling with replacement
                #dstrat<-svydesign(id=~1,strata=~bstrat, data=dat, fpc=~fpc) # stratified sampling without replacement
                #fit=svyglm(as.formula("y~z1+z2+x"%.%i), design=dstrat, family="binomial")
                # option 2: use survey package (two-phase sampling; with and without fpc produce the same results) 
                #dstrat<-twophase(id=list(~1,~1), strata=list(NULL,~bstrat), subset=~ph2, data=dat, fpc=list(NULL,~fpc)) # with fpc
                dstrat<-twophase(id=list(~1,~1), strata=list(NULL,~bstrat), subset=~ph2, data=dat) # without fpc
                fit=svyglm(as.formula("y~z1+z2+x"%.%i), design=dstrat, family="binomial")
                # option 3: use glm without weights
                #fit=glm(as.formula("y~z1+z2+x"%.%i), dat, family=binomial)
            }
            #last(summary(fit)$coef) # for glm
            summary(fit)$coef["x"%.%i,"Pr(>|t|)"] # for svyglm
        })
        p.adj=p.adjust(pvals, method="BH")
        out=min(p.adj)        
      
    } else if (proj=="test1") {
        pvals=sapply (1:p, function(i) {
            fit=glm(as.formula("y~z1+z2+x"%.%i), dat, family=binomial, weights=dat$wt)
            last(summary(fit)$coef)# p 
        })
        out=pvals
      
    } else if (proj=="test2") {
        ests=sapply (1:p, function(i) {
            fit=glm(as.formula("y~z1+z2+x"%.%i), dat, family=binomial, weights=dat$wt)
            last(summary(fit)$coef[,1]) # est
        })
        out=ests
      
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
            
        do.est=function(dat.b){# dat.b has dat.b$case and dat.b$control
        
            use.w= endsWith(proj,"w")
            
            if(startsWith(proj,"perm_min")) {
                dat.train=rbind(data.frame(y=1,dat.b$case), data.frame(y=0,dat.b$control))
                pvals=sapply (1:p, function(i) {
                    set.seed(123)
                    if(!use.w) {
                        # option 3: use glm without weights
                        fit=glm(as.formula("y~z1+z2+x"%.%i), dat.train, family=binomial)
                    } else {
                        # option 1: use survey package
                        #dstrat<-svydesign(id=~1,strata=~bstrat, weights=~wt, data=dat.train) # stratified sampling with replacement
                        #dstrat<-svydesign(id=~1,strata=~bstrat, data=dat.train, fpc=~fpc) # stratified sampling without replacement
                        #fit=svyglm(as.formula("y~z1+z2+x"%.%i), design=dstrat, family="binomial")
                        # option 2: use survey package (two-phase sampling; with and without fpc produce the same results)
                        #dstrat<-twophase(id=list(~1,~1), strata=list(NULL,~bstrat), subset=~ph2, data=dat.train, fpc=list(NULL,~fpc)) # with fpc
                        dstrat<-twophase(id=list(~1,~1), strata=list(NULL,~bstrat), subset=~ph2, data=dat.train) # without fpc
                        fit=svyglm(as.formula("y~z1+z2+x"%.%i), design=dstrat, family="binomial")
                    }
                    #last(summary(fit)$coef) # glm
                    summary(fit)$coef["x"%.%i,"Pr(>|t|)"] # svyglm
                })
                min(pvals)
                
            } else if(startsWith(proj,"perm_rf")) {
            
                splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)                
                
                # no screening
                if(verbose) print("without screening")
                cv.aucs <-  sapply( splits, function(split){
                    dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
                    dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
                    set.seed(123)
                    if(use.w) {
                        # option 2: fit model with weights (weighted)
                        #fit.rf <- ranger( as.formula(f), data = dat.train, case.weights = dat.train$wt, probability = TRUE, min.node.size = 1 )
                        # option 3: fit model without weights (semi-weighted)
                        fit.rf <- ranger( as.formula(f), data=dat.train, probability=TRUE, min.node.size=1 )
                        pred.rf <- predict( fit.rf, data=dat.test )# not allow type="prob"
                        WeightedAUC(WeightedROC(guess=pred.rf$predictions[,'1'], label=dat.test$Y, weight=dat.test$wt))
                    } else {
                        # option 1: fit model with weights (unweighted)
                        fit.rf <- randomForest( as.formula(f), dat.train )# randomForest does not handle weights
                        pred.rf <- predict( fit.rf, newdata=dat.test, type="prob" )
                        fast.auc( pred.rf[,2], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
                    }
                })
                ret=mean(cv.aucs)
                
                if (startsWith(proj,"perm_rf2")) {
                    # ulr screening, p threshold 0.1 (for rv144 cutoff .1 works better than .05)
                    # this learner is about twice as fast as no screening since there are less variables in the learner
                    if(verbose) print("with screening")
                    cv.aucs.2 <-  sapply( splits, function(split){
                        dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
                        dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
                        set.seed(123)
                        
                        # screening
                        predictors=names(dat.train)[startsWith(names(dat.train),"z") | startsWith(names(dat.train),"x")]
                        if(use.w) {
                            sel=screen_ulr (Y=dat.train$Y, X=dat.train[,predictors], family=binomial(), cutoff=0.1, obsWeights=dat.train$wt)
                        } else {
                            sel=screen_ulr (Y=dat.train$Y, X=dat.train[,predictors], family=binomial(), cutoff=0.1)
                        }
                        predictors=predictors[sel]
                        if(verbose) myprint(sum(sel))
                        
                        if (length(predictors)>0) {
                            if(use.w) {
                                fit.rf <- ranger( as.formula("factor(Y)~"%.% concatList(predictors,"+")), data = dat.train, case.weights = dat.train$wt, probability = TRUE, min.node.size = 1 )
                                pred.rf <- predict( fit.rf, data = dat.test )
                                WeightedAUC(WeightedROC(guess=pred.rf$predictions[,'1'], label=dat.test$Y, weight=dat.test$wt))
                            } else {
                                fit.rf <- randomForest( as.formula("factor(Y)~"%.% concatList(predictors,"+")), dat.train )
                                pred.rf <- predict( fit.rf, newdata=dat.test, type="prob" )
                                fast.auc( pred.rf[,2], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
                            }                            
                        } else .5
                    })
                    ret.2=mean(cv.aucs.2)
                    # combining with and without screening by taking maximum
                    ret=max(ret, ret.2)              
                }
                
                if (startsWith(proj,"perm_rf4")) {
                    # add lasso screening, no weighting
                    cv.aucs.2 <-  sapply( splits, function(split){
                        dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
                        dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
                        set.seed(123)
                        # lasso screening
                        predictors=setdiff(names(dat.train),"Y")
                        sel=screen_lasso (Y=dat.train$Y, X=dat.train[-1], family=binomial(), alpha=1) 
                        predictors=predictors[sel]
                        if (length(predictors)>0) {
                            fit.rf <- randomForest( as.formula("factor(Y)~"%.% concatList(predictors,"+")), dat.train )
                            pred.rf <- predict( fit.rf, newdata=dat.test, type="prob" )
                            fast.auc( pred.rf[,2], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
                        } else .5
                    })
                    ret.2=mean(cv.aucs.2)
                    # maximum of with and without screening
                    ret=max(ret, ret.2)              
                }
                
                ret
    
            } else if(startsWith(proj,"perm_dl")) {
                X=as.matrix(rbind(dat.b$case, dat.b$control))
                y=c(rep(1,nrow(dat.b$case)), rep(0,nrow(dat.b$control))) # for run_lgb, y needs to be a vector; for run_dl, y needs to be a matrix
                
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
                X=as.matrix(rbind(dat.b$case, dat.b$control))
                y=c(rep(1,nrow(dat.b$case)), rep(0,nrow(dat.b$control))) # for run_lgb, y needs to be a vector; for run_dl, y needs to be a matrix
                
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
    
            } else if(startsWith(proj,"perm_mlr")) {
                # MLR : Multivariate Logistic Regression #
                splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)                
        
                cv.aucs <-  sapply( splits, function(split){
                  dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
                  dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
                  set.seed(123)
                  if(use.w) {
                    # option 2: use glm without weights (semi-weighted)
                    #fit.mlr <- glm( as.formula(f), dat.train, family=binomial)
                    #pred.mlr <- predict( fit.mlr, newdata=dat.test )
                    # option 3: use survey package (two-phase sampling) (weighted)
                    dstrat<-twophase(id=list(~1,~1), strata=list(NULL,~bstrat), subset=~ph2, data=dat.train)
                    fit.mlr <- svyglm( as.formula(f), design=dstrat, family="binomial")
                    fit.mlr$coefficients[is.na(fit.mlr$coefficients)] <- 0
                    dat.test.temp <- as.matrix(select(dat.test,-c('Y','y','bstrat','ph2','wt','fpc')))
                    pred.mlr <- cbind(1,dat.test.temp)%*%as.matrix(fit.mlr$coefficients)
                    WeightedAUC(WeightedROC(guess=pred.mlr, label=dat.test$Y, weight=dat.test$wt))
                  } else {
                    # option 1: use glm without weights (unweighted)
                    fit.mlr <- glm( as.formula(f), dat.train, family=binomial )
                    pred.mlr <- predict( fit.mlr, newdata=dat.test )
                    fast.auc( pred.mlr, dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
                  }
                })
                mean(cv.aucs)
        
            } else if(startsWith(proj,"perm_rr")) {
                # RR : Ridge Regression #
                splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)                
        
                cv.aucs <-  sapply( splits, function(split){
                  dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
                  dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
                  dat.train.X <- as.matrix(select( dat.train, -c('Y','y','bstrat','ph2','wt','fpc') )); dat.train.y <- as.matrix(select( dat.train, 'Y' ))
                  dat.test.X <- as.matrix(select( dat.test, -c('Y','y','bstrat','ph2','wt','fpc') )); dat.test.y <- as.matrix(select( dat.test, 'Y' ))
                  set.seed(123)
                  if(use.w) {
                    # option 2: fit model without weights (semi-weighted)
                    #fit.rr <- cv.glmnet( x=dat.train.X, y=dat.train.y, family="binomial", type.measure='auc', nfolds=5, alpha=0 ) # Ridge penalty
                    #pred.rr <- as.numeric(drop( predict(fit.rr, newx=dat.test.X, s=fit.rr$lambda.min)))
                    # option 3: fit model with weights (weighted)
                    fit.rr <- cv.glmnet( x=dat.train.X, y=dat.train.y, weights=dat.train$wt, family="binomial", type.measure='auc', nfolds=5, alpha=0 ) # Ridge penalty
                    WeightedAUC(WeightedROC(guess=pred.rr, label=dat.test.y, weight=dat.test$wt))
                  } else {
                    # option 1: fit model without weights (unweighted)
                    fit.rr <- cv.glmnet( x=dat.train.X, y=dat.train.y, family="binomial", type.measure='auc', nfolds=5, alpha=0 ) # Ridge penalty
                    pred.rr <- as.numeric(drop( predict( fit.rr, newx=dat.test.X, s=fit.rr$lambda.min ) ))
                    fast.auc( pred.rr, dat.test.y, reverse.sign.if.nece = FALSE, quiet = TRUE )
                  }
                })
                mean(cv.aucs)
        
            } else if(startsWith(proj,"perm_lr")) {
                # LR : Lasso Regression #
                splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)                
        
                cv.aucs <-  sapply( splits, function(split){
                  dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
                  dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
                  dat.train.X <- as.matrix(select( dat.train, -c('Y','y','bstrat','ph2','wt','fpc') )); dat.train.y <- as.matrix(select(dat.train, 'Y'))
                  dat.test.X <- as.matrix(select( dat.test, -c('Y','y','bstrat','ph2','wt','fpc') )); dat.test.y <- as.matrix(select(dat.test, 'Y'))
                  set.seed(123)
                  if(use.w) {
                    # option 2: fit model without weights (semi-weighted)
                    #fit.lr <- cv.glmnet( x=dat.train.X, y=dat.train.y, family="binomial" , type.measure='auc', nfolds=5, alpha=1 ) # Lasso penalty
                    #pred.lr <- as.numeric(drop( predict( fit.lr , newx = dat.test.X , s = fit.lr$lambda.min ) ))
                    # option 3: fit model with weights (weighted)
                    fit.lr <- cv.glmnet( x=dat.train.X, y=dat.train.y, weights=dat.train$wt, family="binomial", type.measure='auc', nfolds=5, alpha=1 ) # Lasso penalty
                    WeightedAUC(WeightedROC(guess=pred.lr, label=dat.test.y, weight=dat.test$wt))
                  } else {
                    # option 1: fit model without weights (unweighted)
                    fit.lr <- cv.glmnet( x=dat.train.X, y=dat.train.y, family="binomial", type.measure='auc', nfolds=5, alpha=1 ) # Lasso penalty
                    pred.lr <- as.numeric(drop( predict( fit.lr, newx=dat.test.X, s=fit.lr$lambda.min ) ))
                    fast.auc( pred.lr, dat.test.y, reverse.sign.if.nece=FALSE, quiet=TRUE )
                  }
                })
                mean(cv.aucs)
        
            } else if(startsWith(proj,"perm_dt")) {
                # DT : Decision Tree #
                splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)                
        
                cv.aucs <-  sapply( splits, function(split){
                  dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
                  dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
                  set.seed(123)
                  if(use.w) {
                    # option 2: fit model without weights (semi-weighted)
                    #fit.dt <- rpart( as.formula(f), data=dat.train )
                    #pred.dt <- predict( fit.dt, newdata=dat.test, type="prob" )
                    # option 3: fit model with weights (weighted)
                    fit.dt <- rpart( as.formula(f), data=dat.train, weights=dat.train$wt )
                    WeightedAUC(WeightedROC(guess=pred.dt[,2], label=dat.test$Y, weight=dat.test$wt))
                  } else {
                    # option 1: fit model without weights (unweighted)
                    fit.dt <- rpart( as.formula(f), data=dat.train )
                    pred.dt <- predict( fit.dt, newdata=dat.test, type="prob" )
                    fast.auc( pred.dt[,2], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
                  }
                })
                mean(cv.aucs)
        
            } else if(startsWith(proj,"perm_bg")) {
                # BG : Bagging #
                splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)                
        
                cv.aucs <-  sapply( splits, function(split){
                  dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
                  dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
                  set.seed(123)
                  if(use.w) {
                    # option 2: fit model without weights (semi-weighted)
                    #fit.bg <- ranger( as.formula(f), data=dat.train, probability=TRUE, min.node.size=1, mtry=ncol(select(dat.train, -c('Y','y','bstrat','ph2','wt','fpc'))) )
                    #pred.bg <- predict( fit.bg, data = dat.test )# not allow type="prob"
                    # option 3: fit model with weights (weighted)
                    fit.bg <- ranger( as.formula(f), data=dat.train, case.weights=dat.train$wt, probability=TRUE, min.node.size=1, mtry=ncol(select(dat.train, -c('Y','y','bstrat','ph2','wt','fpc'))) )
                    WeightedAUC(WeightedROC(guess=pred.bg$predictions[,'1'], label=dat.test$Y, weight=dat.test$wt))         
                  } else {
                    # option 1: fit model without weights (unweighted)
                    fit.bg <- randomForest( as.formula(f), dat.train, mtry = ncol(select(dat.train, -c('Y','y','bstrat','ph2','wt','fpc'))) )# randomForest does not handle weights
                    pred.bg <- predict( fit.bg, newdata=dat.test, type="prob" )
                    fast.auc( pred.bg[,2], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
                  }
                })
                mean(cv.aucs)
        
            } else if(startsWith(proj,"perm_ab")) {
                # AB : AdaBoost #
                splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)                
        
                cv.aucs <-  sapply( splits, function(split){
                  dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
                  dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
                  dat.train$Y <- factor(dat.train$Y)
                  set.seed(123)
                  if(use.w) {
                    # The weights are used to only calculate CV-AUC (semi-weighted)
                    fit.ab <- adam2( Y~., data=select(dat.train, -c('y','bstrat','ph2','wt','fpc')), size=100, alg='cart' )
                    pred.ab <- predict( fit.ab, newdata=dat.test )
                    WeightedAUC(WeightedROC(guess=pred.ab, label=dat.test$Y, weight=dat.test$wt))
                  } else {
                    # option 1: fit model without weights (unweighted)
                    fit.ab <- adam2( Y~., data=select(dat.train, -c('y','bstrat','ph2','wt','fpc')), size=100, alg='cart' )
                    pred.ab <- predict( fit.ab, newdata=dat.test )
                    fast.auc( pred.ab, dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
                  }
                })
                mean(cv.aucs)
        
            } else if(startsWith(proj,"perm_ubg")) {
                # UBG : Bagging with under-sampling #
                splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)                
        
                cv.aucs <-  sapply( splits, function(split){
                  dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
                  dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
                  class.dist <- table(dat.train$Y)
                  weights.samp <- rep(NA, nrow(dat.train))
                  weights.samp[dat.train$Y==1] <- class.dist['0']/class.dist['1']; weights.samp[dat.train$Y==0] <- 1
                  set.seed(123)
                  if(use.w) {
                    # The weights are used to only calculate CV-AUC (semi-weighted)
                    fit.ubg <- ranger( as.formula(f), data=dat.train, case.weights=weights.samp, probability=TRUE, min.node.size=1, mtry=ncol(select(dat.train, -c('Y','y','bstrat','ph2','wt','fpc'))),
                                       sample.fraction = (class.dist['1']*2)/sum(class.dist) )
                    pred.ubg <- predict( fit.ubg, data=dat.test )
                    WeightedAUC(WeightedROC(guess=pred.ubg$predictions[,'1'], label=dat.test$Y, weight=dat.test$wt))
                  } else {
                    # option 1: fit model without weights (unweighted)
                    fit.ubg <- ranger( as.formula(f), data=dat.train, case.weights=weights.samp, probability=TRUE, min.node.size=1, mtry=ncol(select(dat.train, -c('Y','y','bstrat','ph2','wt','fpc'))),
                                       sample.fraction = (class.dist['1']*2)/sum(class.dist) )
                    pred.ubg <- predict( fit.ubg, data=dat.test )
                    fast.auc( pred.ubg$predictions[,'1'], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
                  }
                })
                mean(cv.aucs)
        
            } else if(startsWith(proj,"perm_obg")) {
                # OBG : Bagging with over-sampling #
                splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)                
        
                cv.aucs <-  sapply( splits, function(split){
                  dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
                  dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
          
                  # Oversampling minority class
                  set.seed( 123 )
                  idx1 <- sample( 1:sum(dat.train$Y==1), sum(dat.train$Y==0), replace = TRUE )
                  dat.train.over <- rbind(dat.train[idx1,,drop=F], subset(dat.train, Y==0))
          
                  set.seed(123)
                  if(use.w) {
                    # The weights are used to only calculate CV-AUC (semi-weighted)
                    fit.obg <- ranger( as.formula(f), data=dat.train.over, probability=TRUE, min.node.size=1, mtry=ncol(select(dat.train.over, -c('Y','y','bstrat','ph2','wt','fpc'))) )
                    pred.obg <- predict( fit.obg, data=dat.test )
                    WeightedAUC(WeightedROC(guess=pred.obg$predictions[,'1'], label=dat.test$Y, weight=dat.test$wt))
                  } else {
                    # option 1: fit model without weights (unweighted)
                    fit.obg <- ranger( as.formula(f), data=dat.train.over, probability=TRUE, min.node.size=1, mtry=ncol(select(dat.train.over, -c('Y','y','bstrat','ph2','wt','fpc'))) )
                    pred.obg <- predict( fit.obg, data=dat.test )
                    fast.auc( pred.obg$predictions[,'1'], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
                  }
                })
                mean(cv.aucs)
        
            } else if(startsWith(proj,"perm_urf")) {
                # URF : RF with under-sampling #
                splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)                
        
                cv.aucs <-  sapply( splits, function(split){
                  dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
                  dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
                  class.dist <- table(dat.train$Y)
                  weights.samp <- rep(NA, nrow(dat.train))
                  weights.samp[dat.train$Y==1] <- class.dist['0']/class.dist['1']; weights.samp[dat.train$Y==0] <- 1
                  set.seed(123)
                  if(use.w) {
                    # The weights are used to only calculate CV-AUC (semi-weighted)
                    fit.urf <- ranger( as.formula(f), data=dat.train, case.weights=weights.samp, probability=TRUE, min.node.size=1, sample.fraction = (class.dist['1']*2)/sum(class.dist) )
                    pred.urf <- predict( fit.urf, data=dat.test )
                    WeightedAUC(WeightedROC(guess=pred.urf$predictions[,'1'], label=dat.test$Y, weight=dat.test$wt))
                  } else {
                    # option 1: fit model without weights (unweighted)
                    fit.urf <- ranger( as.formula(f), data=dat.train, case.weights=weights.samp, probability=TRUE, min.node.size=1, sample.fraction = (class.dist['1']*2)/sum(class.dist) )
                    pred.urf <- predict( fit.urf, data=dat.test )
                    fast.auc( pred.urf$predictions[,'1'], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
                  }
                })
                mean(cv.aucs)
        
            } else if(startsWith(proj,"perm_orf")) {
                # ORF : RF with over-sampling #
                splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)                
        
                cv.aucs <-  sapply( splits, function(split){
                  dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
                  dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
          
                  # Oversampling minority class
                  set.seed( 123 )
                  idx1 <- sample( 1:sum(dat.train$Y==1), sum(dat.train$Y==0), replace = TRUE )
                  dat.train.over <- rbind(dat.train[idx1,,drop=F], subset(dat.train, Y==0))
          
                  set.seed(123)
                  if(use.w) {
                    # The weights are used to only calculate CV-AUC (semi-weighted)
                    fit.orf <- ranger( as.formula(f), data=dat.train.over, probability=TRUE, min.node.size=1 )
                    pred.orf <- predict( fit.orf, data=dat.test )
                    WeightedAUC(WeightedROC(guess=pred.orf$predictions[,'1'], label=dat.test$Y, weight=dat.test$wt))
                  } else {
                    # option 1: fit model without weights (unweighted)
                    fit.orf <- ranger( as.formula(f), data=dat.train.over, probability=TRUE, min.node.size=1 )
                    pred.orf <- predict( fit.orf, data=dat.test )
                    fast.auc( pred.orf$predictions[,'1'], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
                  }
                })
                mean(cv.aucs)
        
            } else if(startsWith(proj,"perm_rusb")) {
                # RUSB : Random Under Sampling Boosting #
                splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)                
        
                cv.aucs <-  sapply( splits, function(split){
                  dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
                  dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
                  dat.train$Y <- factor(dat.train$Y)
                  set.seed(123)
                  if(use.w) {
                    # The weights are used to only calculate CV-AUC (semi-weighted)
                    fit.rusb <- rus( Y~., data=select(dat.train, -c('y','bstrat','ph2','wt','fpc')), size=100, alg='cart' )
                    pred.rusb <- predict( fit.rusb, newdata = dat.test )
                    WeightedAUC(WeightedROC(guess=pred.rusb, label=dat.test$Y, weight=dat.test$wt))
                  } else {
                    # option 1: fit model without weights (unweighted)
                    fit.rusb <- rus( Y~., data=select(dat.train, -c('y','bstrat','ph2','wt','fpc')), size=100, alg='cart' )
                    pred.rusb <- predict( fit.rusb, newdata = dat.test )
                    fast.auc( pred.rusb, dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
                  }
                })
                mean(cv.aucs)
        
            } else if(startsWith(proj,"perm_smoteb")) {
                # SMOTEB : SMOTEBoost #
                splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)                
        
                cv.aucs <-  sapply( splits, function(split){
                  dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
                  dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
                  dat.train$Y <- factor(dat.train$Y)
                  set.seed(123)
                  if(use.w) {
                    # The weights are used to only calculate CV-AUC (semi-weighted)
                    fit.smoteb <- sbo( Y~., data=select(dat.train, -c('y','bstrat','ph2','wt','fpc')), size=100, alg='cart', over=500 )
                    pred.smoteb <- predict( fit.smoteb, newdata = dat.test )
                    WeightedAUC(WeightedROC(guess=pred.smoteb, label=dat.test$Y, weight=dat.test$wt))
                  } else {
                    # option 1: fit model without weights (unweighted)
                    fit.smoteb <- sbo( Y~., data=select(dat.train, -c('y','bstrat','ph2','wt','fpc')), size=100, alg='cart', over=500 )
                    pred.smoteb <- predict( fit.smoteb, newdata = dat.test )
                    fast.auc( pred.smoteb, dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
                  }
                })
                mean(cv.aucs)
        
            } else stop("wrong proj")
        }

        dat.b=list(case=subset(dat, y==1), control=subset(dat, y==0))
        est=do.est(dat.b)                
        
        dat.case.control=rbind(dat.b$case, dat.b$control) # for permuting 
        if (choose(N,n1)<=nperm) p.method="exact" else p.method="Monte Carlo"
        if(p.method=="exact"){        
            indexes=combn(1:N, n1) 
            ref.distr=apply(indexes, 2, function(ind.1) {
                ind.0=setdiff(1:N, ind.1)
                dat.b=list()
                dat.b$case=   cbind(dat.case.control[1:n1,   c("z1","z2","bstrat","wt","ph2","fpc")], subset(dat.case.control[ind.1,], select=-c(z1,z2,bstrat,wt,ph2,fpc)) )
                dat.b$control=cbind(dat.case.control[1:n0+n1,c("z1","z2","bstrat","wt","ph2","fpc")], subset(dat.case.control[ind.0,], select=-c(z1,z2,bstrat,wt,ph2,fpc)) )
#                dat.b$case=   dat.case.control[ind.1,]
#                dat.b$control=dat.case.control[ind.0,]
                do.est(dat.b)
            })
            
        } else if (p.method=="Monte Carlo") {
            #### save rng state before set.seed in order to restore before exiting this function
            save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
            if (class(save.seed)=="try-error") { set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }                        
            
            ref.distr=sapply(1:nperm, function(perm.i) {
                if(verbose) print(perm.i)
                set.seed(perm.i)
                ind.1=sample(1:N, n1)    # cases in the permutated dataset
                ind.0=setdiff(1:N,ind.1) # controls in the permutated dataset
                dat.b=list()
                dat.b$case=   cbind(dat.case.control[1:n1,   c("z1","z2","bstrat","wt","ph2","fpc")], subset(dat.case.control[ind.1,], select=-c(z1,z2,bstrat,wt,ph2,fpc)) )
                dat.b$control=cbind(dat.case.control[1:n0+n1,c("z1","z2","bstrat","wt","ph2","fpc")], subset(dat.case.control[ind.0,], select=-c(z1,z2,bstrat,wt,ph2,fpc)) )
                do.est(dat.b)
            })            
                
            #### restore rng state 
            assign(".Random.seed", save.seed, .GlobalEnv)     
        } else stop("wrong p.method")
        
        # plot permutation samples distribution
        if(make.plot) {plot(density(ref.distr), main=seed); abline(v=est, lty=2); abline(v=0.5)}
        
        # get pvalue
        if(startsWith(proj,"perm_min")) {
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


library(parallel); library(kyotil)
source("../helper_cvauc.R")
mc.cores=5
beta=0
res=mclapply(1:1e3, mc.cores=mc.cores, function (seed) {
    dat=sim.rv144(n1=NULL, n0=NULL, seed, alpha=if(beta==1) -7.6 else -6, betas=if(beta==1) c(0.5,0.5, 0.8,0,-2)  else c(0,0,0,0,0,0), beta.z.1=0.5, beta.z.2=0, n=1e4) 
    c(nrow(dat), sum(dat$y))
#    dat.1=sim.rv144(n1=40, n0=200, seed, alpha=-1.9,    betas=c(0.9,0.7,1.2,-1,-1.2), beta.z.1=0.5, beta.z.2=0)        
#    dat.2=sim.rv144(n1=NULL, n0=NULL, seed, alpha=-7.6, betas=c(0.9,0.7,1.2,-1,-1.2), beta.z.1=0.5, beta.z.2=0, n=1e4)   
#    c(nrow(dat.1), sum(dat.1$y), nrow(dat.2), sum(dat.2$y))
})
tmp=do.call(rbind, res)
apply(tmp, 2, mean)
apply(tmp, 2, median)


# rv144 median sample size and case number
# cohort
# 240  50 
# twophase
# 186  31
