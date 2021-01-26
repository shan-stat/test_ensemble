### Import helper function ###
setwd("set your local directory where helper_cvauc.R is located")
source("helper_cvauc.R")



### 1. Setting process arguments ###
Args <- commandArgs(trailingOnly=TRUE) 
if (length(Args)==0) {
  # batch.size and batch.number are need to a batch script
  # sim.setting: rv144ph2_n10000_1 for RV144 experiments in Section 3
  # fit.setting: 5fold for 5-fold cross validation
  # proj: you can choose one of tests listed in Section 3
  # ipw: uw for unweighted; sw for semi-weights; w for weighted
  Args=c(batch.size="10",batch.number="1",sim.setting="rv144ph2_n10000_1",fit.setting="5fold",proj="BH",ipw="uw")
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
i=i+1; proj=Args[i]
i=i+1; ipw=Args[i]

nperm=1e3 # permutation replicates
verbose=ifelse(unix(),0,2)
myprint(sample.size, sim.model)



### 2. Experiments ###
res=sapply(seeds, simplify="array", function (seed) {
  
  myprint(seed)
  
  # Sample size #
  if (sample.size=="n10000") {
    n1=NULL; n0=NULL; n=1e4
  } else stop("wrong simple.size")
  
  # Simulated dataset #
  if (sim.model=="rv144ph2")  {
    if((proj=='BH' & ipw=='w')|(proj=='perm_MP' & ipw=='w')){
      dat=sim.rv144(n1, n0, seed, alpha=if(beta==1) -7.6 else -6, betas=c(0.5,0.5, 0.8,0,-2), beta.z.1=0.5, beta.z.2=0, n=n, full=TRUE) 
    } else {
      dat=sim.rv144(n1, n0, seed, alpha=if(beta==1) -7.6 else -6, betas=c(0.5,0.5, 0.8,0,-2), beta.z.1=0.5, beta.z.2=0, n=n) 
    }
    f = "as.factor(Y)~z1+z2+" %.% concatList(names(dat)[startsWith(names(dat),"x")],"+")    
  } else stop ("wrong sim.model")
  
  # Actual number of cases and controls #
  N=nrow(dat)
  n1=sum(dat$y)
  n0=N-n1
  p=sum(startsWith(names(dat),"x"))
  myprint(n1,n0,p)
  
  # 2-1. Fitting tests #
  if (startsWith(proj,"BH")) {
    # BH: Benjamini-Hochberg #
    pvals=sapply (1:p, function(i) {
      set.seed(123)
      if(ipw=='uw'){
        # Unweighted # 
        fit=glm(as.formula("y~z1+z2+x"%.%i), dat, family=binomial)
        last(summary(fit)$coef)
      } else if(ipw=='w'){
        # Weighted #
        dstrat<-twophase(id=list(~1,~1), strata=list(NULL,~bstrat), subset=~ph2, data=dat)
        fit=svyglm(as.formula("y~z1+z2+x"%.%i), design=dstrat, family="binomial")
        summary(fit)$coef["x"%.%i,"Pr(>|t|)"]
      }
    })
    p.adj=p.adjust(pvals, method="BH")
    out=min(p.adj)
    
  } else if (startsWith(proj,"perm")) {   
    # perm: Permutation-based tests #
    do.est=function(dat.b){
      
      if(startsWith(proj,"perm_MP")) {
        # MP: the minimum p-value #
        dat.train=rbind(data.frame(y=1,dat.b$case), data.frame(y=0,dat.b$control))
        pvals=sapply (1:p, function(i) {
          set.seed(123)
          if(ipw=='uw') {
            # Unweighted #
            fit=glm(as.formula("y~z1+z2+x"%.%i), dat.train, family=binomial)
            last(summary(fit)$coef)
          } else if(ipw=='w') {
            # Weighted #
            dstrat<-twophase(id=list(~1,~1), strata=list(NULL,~bstrat), subset=~ph2, data=dat.train)
            fit=svyglm(as.formula("y~z1+z2+x"%.%i), design=dstrat, family="binomial")
            summary(fit)$coef["x"%.%i,"Pr(>|t|)"]
          }
        })
        min(pvals)
        
      } else if(startsWith(proj,"perm_RR")) {
        # RR: ridge regression #
        splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)                
        
        cv.aucs <-  sapply( splits, function(split){
          dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
          dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
          dat.train.X <- as.matrix(select( dat.train, -c('Y','y','bstrat','ph2','wt','fpc') )); dat.train.y <- as.matrix(select( dat.train, 'Y' ))
          dat.test.X <- as.matrix(select( dat.test, -c('Y','y','bstrat','ph2','wt','fpc') )); dat.test.y <- as.matrix(select( dat.test, 'Y' ))
          set.seed(123)
          if(ipw=='uw'|ipw=='sw') {
            fit.rr <- cv.glmnet( x=dat.train.X, y=dat.train.y, family="binomial", type.measure='auc', nfolds=5, alpha=0 ) # Ridge penalty
            pred.rr <- as.numeric(drop( predict( fit.rr, newx=dat.test.X, s=fit.rr$lambda.min ) ))
            if(ipw=='uw'){
              # Unweighted #
              fast.auc( pred.rr, dat.test.y, reverse.sign.if.nece = FALSE, quiet = TRUE )
            } else if(ipw=='sw'){
              # semi-weighted #
              WeightedAUC(WeightedROC(guess=pred.rr, label=dat.test.y, weight=dat.test$wt))
            }
          } else if(ipw=='w'){
            # Weighted #
            fit.rr <- cv.glmnet( x=dat.train.X, y=dat.train.y, weights=dat.train$wt, family="binomial", type.measure='auc', nfolds=5, alpha=0 ) # Ridge penalty
            pred.rr <- as.numeric(drop( predict(fit.rr, newx=dat.test.X, s=fit.rr$lambda.min)))
            WeightedAUC(WeightedROC(guess=pred.rr, label=dat.test.y, weight=dat.test$wt))
          }
        })
        mean(cv.aucs)
        
      } else if(startsWith(proj,"perm_BG")) {
        # BG: bagging #
        splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)                
        
        cv.aucs <-  sapply( splits, function(split){
          dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
          dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
          set.seed(123)
          if(ipw=='uw'|ipw=='sw') {
            fit.bg <- ranger( as.formula(f), data=dat.train, probability=TRUE, min.node.size=1, mtry=ncol(select(dat.train, -c('Y','y','bstrat','ph2','wt','fpc'))) )
            pred.bg <- predict( fit.bg, data = dat.test )
            if(ipw=='uw'){
              # Unweighted #
              fast.auc( pred.bg$predictions[,'1'], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
            } else if(ipw=='sw'){
              # Semi-weighted #
              WeightedAUC(WeightedROC(guess=pred.bg$predictions[,'1'], label=dat.test$Y, weight=dat.test$wt))
            }
          } else if(ipw=='w'){
            # Weighted #
            fit.bg <- ranger( as.formula(f), data=dat.train, case.weights=dat.train$wt, probability=TRUE, min.node.size=1, mtry=ncol(select(dat.train, -c('Y','y','bstrat','ph2','wt','fpc'))) )
            pred.bg <- predict( fit.bg, data = dat.test )
            WeightedAUC(WeightedROC(guess=pred.bg$predictions[,'1'], label=dat.test$Y, weight=dat.test$wt))
          }
        })
        mean(cv.aucs)
        
      } else if(startsWith(proj,"perm_RF")) {
        # RF: random forest #
        splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)                
        
        cv.aucs <-  sapply( splits, function(split){
          dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
          dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
          set.seed(123)
          if(ipw=='uw'|ipw=='sw') {
            fit.rf <- ranger( as.formula(f), data=dat.train, probability=TRUE, min.node.size=1 )
            pred.rf <- predict( fit.rf, data=dat.test )
            if(ipw=='uw'){
              # Unweighted #
              fast.auc( pred.rf$predictions[,'1'], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
            } else if(ipw=='sw'){
              # Semi-weighted #
              WeightedAUC(WeightedROC(guess=pred.rf$predictions[,'1'], label=dat.test$Y, weight=dat.test$wt))
            }
          } else if(ipw=='w'){
            # Weighted #
            fit.rf <- ranger( as.formula(f), data = dat.train, case.weights = dat.train$wt, probability = TRUE, min.node.size = 1 )
            pred.rf <- predict( fit.rf, data=dat.test )
            WeightedAUC(WeightedROC(guess=pred.rf$predictions[,'1'], label=dat.test$Y, weight=dat.test$wt))
          }
        })
        ret=mean(cv.aucs)
        ret
        
      } else if(startsWith(proj,"perm_AB")) {
        # AB: adaptive boosting #
        splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)                
        
        cv.aucs <-  sapply( splits, function(split){
          dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
          dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
          dat.train$Y <- factor(dat.train$Y)
          set.seed(123)
          fit.ab <- adam2( Y~., data=select(dat.train, -c('y','bstrat','ph2','wt','fpc')), size=100, alg='cart' )
          pred.ab <- predict( fit.ab, newdata=dat.test )
          if(ipw=='uw') {
            # Unweighted #
            fast.auc( pred.ab, dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
          } else if(ipw=='sw'){
            # Semi-weighted #
            WeightedAUC(WeightedROC(guess=pred.ab, label=dat.test$Y, weight=dat.test$wt))
          }
        })
        mean(cv.aucs)
        
      } else if(startsWith(proj,"perm_underBG")) {
        # underBG: bagging with under-sampling #
        splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)                
        
        cv.aucs <-  sapply( splits, function(split){
          dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
          dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
          class.dist <- table(dat.train$Y)
          weights.samp <- rep(NA, nrow(dat.train))
          weights.samp[dat.train$Y==1] <- class.dist['0']/class.dist['1']; weights.samp[dat.train$Y==0] <- 1
          set.seed(123)
          fit.ubg <- ranger( as.formula(f), data=dat.train, case.weights=weights.samp, probability=TRUE, min.node.size=1, mtry=ncol(select(dat.train, -c('Y','y','bstrat','ph2','wt','fpc'))),
                             sample.fraction = (class.dist['1']*2)/sum(class.dist) )
          pred.ubg <- predict( fit.ubg, data=dat.test )
          if(ipw=='uw') {
            # Unweighted #
            fast.auc( pred.ubg$predictions[,'1'], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
          } else if(ipw=='sw'){
            # Semi-weighted #
            WeightedAUC(WeightedROC(guess=pred.ubg$predictions[,'1'], label=dat.test$Y, weight=dat.test$wt))
          }
        })
        mean(cv.aucs)
        
      } else if(startsWith(proj,"perm_underRF")) {
        # underRF: random forest with under-sampling #
        splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)                
        
        cv.aucs <-  sapply( splits, function(split){
          dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
          dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
          class.dist <- table(dat.train$Y)
          weights.samp <- rep(NA, nrow(dat.train))
          weights.samp[dat.train$Y==1] <- class.dist['0']/class.dist['1']; weights.samp[dat.train$Y==0] <- 1
          set.seed(123)
          fit.urf <- ranger( as.formula(f), data=dat.train, case.weights=weights.samp, probability=TRUE, min.node.size=1, sample.fraction = (class.dist['1']*2)/sum(class.dist) )
          pred.urf <- predict( fit.urf, data=dat.test )
          if(ipw=='uw') {
            # Unweighted #
            fast.auc( pred.urf$predictions[,'1'], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
          } else if(ipw=='sw'){
            # Semi-weighted #
            WeightedAUC(WeightedROC(guess=pred.urf$predictions[,'1'], label=dat.test$Y, weight=dat.test$wt))
          }
        })
        mean(cv.aucs)
        
      } else if(startsWith(proj,"perm_underAB")) {
        # underAB: adaboost with under-sampling #
        splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)                
        
        cv.aucs <-  sapply( splits, function(split){
          dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
          dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
          dat.train$Y <- factor(dat.train$Y)
          set.seed(123)
          fit.rusb <- rus( Y~., data=select(dat.train, -c('y','bstrat','ph2','wt','fpc')), size=100, alg='cart' )
          pred.rusb <- predict( fit.rusb, newdata = dat.test )
          if(ipw=='uw') {
            # Unweighted #
            fast.auc( pred.rusb, dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
          } else if(ipw=='sw'){
            # Semi-weighted #
            WeightedAUC(WeightedROC(guess=pred.rusb, label=dat.test$Y, weight=dat.test$wt))
          }
        })
        mean(cv.aucs)
        
      } else if(startsWith(proj,"perm_overBG")) {
        # overBG: bagging with over-sampling #
        splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)                
        
        cv.aucs <-  sapply( splits, function(split){
          dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
          dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
          # over-sampling minority class #
          set.seed(123)
          idx1 <- sample( 1:sum(dat.train$Y==1), sum(dat.train$Y==0), replace = TRUE )
          dat.train.over <- rbind(dat.train[idx1,,drop=F], subset(dat.train, Y==0))
          set.seed(123)
          fit.obg <- ranger( as.formula(f), data=dat.train.over, probability=TRUE, min.node.size=1, mtry=ncol(select(dat.train.over, -c('Y','y','bstrat','ph2','wt','fpc'))) )
          pred.obg <- predict( fit.obg, data=dat.test )
          if(ipw=='uw') {
            # Unweighted #
            fast.auc( pred.obg$predictions[,'1'], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
          } else if(ipw=='sw'){
            # Semi-weighted #
            WeightedAUC(WeightedROC(guess=pred.obg$predictions[,'1'], label=dat.test$Y, weight=dat.test$wt))
          }
        })
        mean(cv.aucs)
        
      } else if(startsWith(proj,"perm_overRF")) {
        # underRF: random forest with over-sampling #
        splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)                
        
        cv.aucs <-  sapply( splits, function(split){
          dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
          dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
          # over-sampling minority class #
          set.seed(123)
          idx1 <- sample( 1:sum(dat.train$Y==1), sum(dat.train$Y==0), replace = TRUE )
          dat.train.over <- rbind(dat.train[idx1,,drop=F], subset(dat.train, Y==0))
          set.seed(123)
          fit.orf <- ranger( as.formula(f), data=dat.train.over, probability=TRUE, min.node.size=1 )
          pred.orf <- predict( fit.orf, data=dat.test )
          if(ipw=='uw') {
            # Unweighted #
            fast.auc( pred.orf$predictions[,'1'], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
          } else if(ipw=='sw'){
            # Semi-weighted #
            WeightedAUC(WeightedROC(guess=pred.orf$predictions[,'1'], label=dat.test$Y, weight=dat.test$wt))
          }
        })
        mean(cv.aucs)
        
      } else if(startsWith(proj,"perm_overAB")) {
        # overAB: adaboost with over-sampling #
        splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)                
        
        cv.aucs <-  sapply( splits, function(split){
          dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
          dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
          dat.train$Y <- factor(dat.train$Y)
          set.seed(123)
          fit.smoteb <- sbo( Y~., data=select(dat.train, -c('y','bstrat','ph2','wt','fpc')), size=100, alg='cart', over=500 )
          pred.smoteb <- predict( fit.smoteb, newdata = dat.test )
          if(ipw=='uw') {
            # Unweighted #
            fast.auc( pred.smoteb, dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
          } else if(ipw=='sw'){
            # Semi-weighted #
            WeightedAUC(WeightedROC(guess=pred.smoteb, label=dat.test$Y, weight=dat.test$wt))
          }
        })
        mean(cv.aucs)
        
      } else stop("wrong proj")
    }
    
    # Estimated CV-AUC #
    dat.b=list(case=subset(dat, y==1), control=subset(dat, y==0))
    est=do.est(dat.b)                
    
    # Permutation process #
    dat.case.control=rbind(dat.b$case, dat.b$control) # for permuting 
    # Save rng state before set.seed in order to restore before exiting this function #
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
    
    # Restore rng state #
    assign(".Random.seed", save.seed, .GlobalEnv)     
    
    # P-value #
    if(startsWith(proj,"perm_MP")) {
      p.value=mean(ref.distr<est)
    } else {
      p.value=mean(ref.distr>=est)
    } 
    
    out=c(est=est, p.value=p.value)
  } else stop("wrong proj")
  
  out
})
res



### 3. Figure B.1 in supplementary material ###
dat=sim.rv144(n1=NULL, n0=NULL, seed=1, alpha=-7.6, betas=c(0.5,0.5, 0.8,0,-2), beta.z.1=0.5, beta.z.2=0, n=10000)
cor=cor(select(dat, -c(y,z1,z2,bstrat,ph2,wt,fpc)))
breaks=c(-1,-.7,-.5,-.3,-.1,.1,.3,.5,.7,1)
hU=DMHeatMap(cor, trace="none", symm=T, dendrogram="none", col=RColorBrewer::brewer.pal(length(breaks)-1,"RdYlGn"), 
             distfun = function(c) as.dist(1 - c), axis=F, breaks=breaks, margins=c(0,0), 
             key = F, Rowv=NA, lower.left.only=FALSE, heatmapOnly=T) 
