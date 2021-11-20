##### Section 2.2, 2.3, and Supplementary Materials Section D #####
source("helper_cvauc.R")
# the following scripts are modified from the caretEnsembleR package
source('caretEnsembleR/helper_functions.R')
source('caretEnsembleR/caretList.R')
source('caretEnsembleR/caretEnsemble.R')
source('caretEnsembleR/caretStack.R')


### 1. Setting arguments ###
Args <- commandArgs(trailingOnly=TRUE) 
if (length(Args)==0) {
  # batch.size and batch.number are for running a batch script
  # sim.setting: "sim.rv144" for simulation study with normal distribution in Section 2.2 and 2.3
  #            : "sim.rv144.lognorm" for simulation study with log-normal distribution in Supplementary Materials Section D
  #            : _1 for power study and _0 for size study
  # sim.linear: "TRUE" for linear scenario; "FALSE" for nonlinear scenario
  # fit.setting: "5-fold" by default
  # proj: "perm_WY" (WY), "perm_RR", "perm_BG", "perm_RF", "perm_AB", "perm_ST" (perm_stacking)
  # ipw: "uw" for unweighted; "sw" for semi-weighted; "w" for weighted
  # scr: "TRUE" for lasso screening; "FALSE" for no screening
  Args=c(batch.size="1000",batch.number="1",sim.setting="sim.rv144_0",sim.linear="FALSE",fit.setting="5fold",proj="perm_WY",ipw="uw",scr="FALSE")
}
myprint(Args)
i=0;
i=i+1; batch.size=as.numeric(Args[i])
i=i+1; batch=as.numeric(Args[i])
seeds=1:batch.size+batch.size*(batch-1); names(seeds)=seeds
myprint(batch, batch.size)
i=i+1; sim.setting=Args[i]; tmp=strsplit(sim.setting,"_")[[1]]
sim.model=tmp[1]
beta=as.numeric(tmp[2]) 
i=i+1; sim.linear=Args[i]
i=i+1; fit.setting=Args[i]; cv.scheme=fit.setting
i=i+1; proj=Args[i]
i=i+1; ipw=Args[i]
i=i+1; scr=Args[i]

nperm=1e3 # permutation replicates
verbose=ifelse(unix(),0,2)

### 2. Experiments ###
res=sapply(seeds, simplify="array", function (seed) {
  
  myprint(seed)
  
  # Sample size #
  n1=NULL; n0=NULL; n=1e4
  
  # Simulated dataset
  if (sim.model=="sim.rv144")  {
      if(sim.linear=="TRUE"){
        dat=sim.rv144(n1, n0, seed, alpha=if(beta==1) -7.4 else -6.1, betas=if(beta==1) c(1.2,1.2,1,-0.8,0) else c(0,0,0,0,0), beta.z.1=0.5, beta.z.2=0, n=n) # for linear
      } else {
        dat=sim.rv144(n1, n0, seed, alpha=if(beta==1) -7.5 else -6.1, betas=if(beta==1) c(0.6,0.6,0.4,0,-2) else c(0,0,0,0,0), beta.z.1=0.5, beta.z.2=0, n=n) # for nonlinear
      }
  } else if (sim.model=="sim.rv144.lognorm")  {
      if(sim.linear=="TRUE"){
        dat=sim.rv144.lognorm(n1, n0, seed, alpha=if(beta==1) -7.4 else -6.1, betas=if(beta==1) c(1.2,1.2,1,-0.8,0) else c(0,0,0,0,0), beta.z.1=0.5, beta.z.2=0, n=n) # for linear
      } else {
        dat=sim.rv144.lognorm(n1, n0, seed, alpha=if(beta==1) -7.5 else -6.1, betas=if(beta==1) c(0.6,0.6,0.4,0,-2) else c(0,0,0,0,0), beta.z.1=0.5, beta.z.2=0, n=n) # for nonlinear
      }
  } else stop ("wrong sim.model")
  
  # Actual number of cases and controls #
  N=nrow(dat)
  n1=sum(dat$y)
  n0=N-n1
  p=sum(startsWith(names(dat),"x"))
  myprint(n1,n0,p)
  f = "as.factor(Y)~z1+z2+" %.% concatList(names(dat)[startsWith(names(dat),"x")],"+")
  
  # 2-1. Fitting testing methods  #
  if (startsWith(proj,"perm")) {   
    # perm: Permutation-based tests #
    do.est=function(dat.b){
      
      if(startsWith(proj,"perm_WY")) {
        # WY #
        dat.train=rbind(data.frame(y=1,dat.b$case), data.frame(y=0,dat.b$control))
        pvals=sapply (1:p, function(i) {
          set.seed(123)
          # Unweighted #
          fit=glm(as.formula("y~z1+z2+x"%.%i), dat.train, family=binomial)
          last(summary(fit)$coef)
        })
        min(pvals)
        
      } else if(startsWith(proj,"perm_RR")) {
        # RR: ridge regression #
        splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)
        cv.aucs <-  sapply( splits, function(split){
          dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
          dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
          
          # Screening #
          screened.var <- rep(TRUE,sum(startsWith(names(dat.train),"x"))); names(screened.var) <- names(dat.train)[startsWith(names(dat.train),"x")]
          if(scr=="TRUE"){
            vars <- names(dat.train)[startsWith(names(dat.train),"x")]
            screened.var <- screen_lasso(Y=dat.train$Y, X=dat.train[,vars], family="binomial") # lasso
          }
          
          # Fitting #
          if(scr=="TRUE" & sum(screened.var)==0){
            0.5
          } else {
            dat.train <- select(dat.train,c("Y","y","z1","z2","bstrat",names(which(screened.var)),"wt","ph2","fpc"))
            dat.test <- select(dat.test,c("Y","y","z1","z2","bstrat",names(which(screened.var)),"wt","ph2","fpc"))
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
          }
        })
        mean(cv.aucs)
        
      } else if(startsWith(proj,"perm_BG")) {
        # BG: bagging #
        splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)
        cv.aucs <-  sapply( splits, function(split){
          dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
          dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
          
          # Screening #
          screened.var <- rep(TRUE,sum(startsWith(names(dat.train),"x"))); names(screened.var) <- names(dat.train)[startsWith(names(dat.train),"x")]
          if(scr=="TRUE"){
            vars <- names(dat.train)[startsWith(names(dat.train),"x")]
            screened.var <- screen_lasso(Y=dat.train$Y, X=dat.train[,vars], family="binomial") # lasso
          }
          
          # Fitting #
          if(scr=="TRUE" & sum(screened.var)==0){
            0.5
          } else {
            f = "as.factor(Y)~z1+z2+" %.% concatList(names(which(screened.var)),"+")
            set.seed(123)
            if(ipw=='uw'|ipw=='sw') {
              fit.bg <- randomForest(as.formula(f), data=dat.train, mtry=sum(screened.var)+2)
              pred.bg <- predict(fit.bg, newdata=dat.test, type='prob')
              if(ipw=='uw'){
                # Unweighted #
                fast.auc( pred.bg[,2], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
              } else if(ipw=='sw'){
                # Semi-weighted #
                WeightedAUC(WeightedROC(guess=pred.bg[,2], label=dat.test$Y, weight=dat.test$wt))
              }
            } else if(ipw=='w'){
              # Weighted #
              fit.bg <- ranger( as.formula(f), data=dat.train, case.weights=dat.train$wt, probability=TRUE, min.node.size=1, mtry=sum(screened.var)+2 )
              pred.bg <- predict( fit.bg, data = dat.test )
              WeightedAUC(WeightedROC(guess=pred.bg$predictions[,'1'], label=dat.test$Y, weight=dat.test$wt))
            }
          }
        })
        mean(cv.aucs)

      } else if(startsWith(proj,"perm_RF")) {
        # RF: random forest #
        splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)
        cv.aucs <-  sapply( splits, function(split){
          dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
          dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
          
          # Screening #
          screened.var <- rep(TRUE,sum(startsWith(names(dat.train),"x"))); names(screened.var) <- names(dat.train)[startsWith(names(dat.train),"x")]
          if(scr=="TRUE"){
            vars <- names(dat.train)[startsWith(names(dat.train),"x")]
            screened.var <- screen_lasso(Y=dat.train$Y, X=dat.train[,vars], family="binomial") # lasso
          }
          
          # Fitting #
          if(scr=="TRUE" & sum(screened.var)==0){
            0.5
          } else {
            f = "as.factor(Y)~z1+z2+" %.% concatList(names(which(screened.var)),"+")
            set.seed(123)
            if(ipw=='uw'|ipw=='sw') {
              fit.rf <- randomForest(as.formula(f), data=dat.train)
              pred.rf <- predict(fit.rf, newdata=dat.test, type="prob")
              if(ipw=='uw'){
                # Unweighted #
                fast.auc( pred.rf[,2], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
              } else if(ipw=='sw'){
                # Semi-weighted #
                WeightedAUC(WeightedROC(guess=pred.rf[,2], label=dat.test$Y, weight=dat.test$wt))
              }
            } else if(ipw=='w'){
              # Weighted #
              fit.rf <- ranger(as.formula(f), data=dat.train, case.weights=dat.train$wt, probability=TRUE, min.node.size=1)
              pred.rf <- predict( fit.rf, data=dat.test )
              WeightedAUC(WeightedROC(guess=pred.rf$predictions[,'1'], label=dat.test$Y, weight=dat.test$wt))
            }
          }
        })
        mean(cv.aucs)
        
      } else if(startsWith(proj,"perm_AB")) {
        # AB: adaptive boosting #
        splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)
        cv.aucs <-  sapply( splits, function(split){
          dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
          dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
          dat.train$Y <- factor(dat.train$Y)
          
          # Screening #
          screened.var <- rep(TRUE,sum(startsWith(names(dat.train),"x"))); names(screened.var) <- names(dat.train)[startsWith(names(dat.train),"x")]
          if(scr=="TRUE"){
            vars <- names(dat.train)[startsWith(names(dat.train),"x")]
            screened.var <- screen_lasso(Y=dat.train$Y, X=dat.train[,vars], family="binomial") # lasso
          }
          
          # Fitting #
          if(scr=="TRUE" & sum(screened.var)==0){
            0.5
          } else {
            set.seed(123)
            fit.ab <- adam2( Y~., data=select(dat.train, c("Y","z1","z2",names(which(screened.var)))), size=100, alg='cart' )
            pred.ab <- predict( fit.ab, newdata=dat.test )
            if(ipw=='uw') {
              # Unweighted #
              fast.auc( pred.ab, dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
            } else if(ipw=='sw'){
              # Semi-weighted #
              WeightedAUC(WeightedROC(guess=pred.ab, label=dat.test$Y, weight=dat.test$wt))
            }
          }
        })
        mean(cv.aucs)
        
      } else if(startsWith(proj,"perm_ST")) {
        # ST: stacking WY and random forest #
        # Selecting a biomarker that has the smallest p value #
        dat.train=rbind(data.frame(y=1,dat.b$case), data.frame(y=0,dat.b$control))
        pvals=sapply (1:p, function(i) {
          set.seed(123)
          fit=glm(as.formula("y~z1+z2+x"%.%i), dat.train, family=binomial)
          last(summary(fit)$coef)
        })
        mp.index <- which.min(pvals)
        
        # For stacking #
        all.var <- rep(TRUE,sum(startsWith(names(dat.train),"x")))
        var.index.set <- list(rf=c(TRUE,TRUE,TRUE,all.var)) # the first three TRUE are for Y, z1, and z2
        set.seed(123)
        pred.st <- get.st.auc(dat.train=dat.train, mp.index=mp.index, var.index=var.index.set, method='method.NNloglik')
        fast.auc( pred.st, dat.train$y, reverse.sign.if.nece = FALSE, quiet = TRUE )
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
      ind.1=sample(1:N, n1)    # cases in the permuted dataset
      ind.0=setdiff(1:N,ind.1) # controls in the permuted dataset
      dat.b=list()
      dat.b$case=   cbind(dat.case.control[1:n1,   c("z1","z2","bstrat","wt","ph2","fpc")], subset(dat.case.control[ind.1,], select=-c(z1,z2,bstrat,wt,ph2,fpc)) )
      dat.b$control=cbind(dat.case.control[1:n0+n1,c("z1","z2","bstrat","wt","ph2","fpc")], subset(dat.case.control[ind.0,], select=-c(z1,z2,bstrat,wt,ph2,fpc)) )
      do.est(dat.b)
    })            
    
    # Restore rng state #
    assign(".Random.seed", save.seed, .GlobalEnv)     
    
    # P-value #
    if(startsWith(proj,"perm_WY")) {
      p.value=mean(ref.distr<est)
    } else {
      p.value=mean(ref.distr>=est)
    } 
    
    out=c(est=est, p.value=p.value)
  } else stop("wrong proj")
  
  out
})
res

### 3. Figure A.2 in the Supporting Information ###
dat=sim.rv144(n1=NULL, n0=NULL, seed=1, alpha=-7.5, betas=c(0.6,0.6,0.4,0,-2), beta.z.1=0.5, beta.z.2=0, n=10000) # for nonlinear
cor=cor(select(dat, -c(y,z1,z2,bstrat,ph2,wt,fpc)))
breaks=c(-1,-.7,-.5,-.3,-.1,.1,.3,.5,.7,1)
hU=DMHeatMap(cor, trace="none", symm=T, dendrogram="none", col=RColorBrewer::brewer.pal(length(breaks)-1,"RdYlGn"), 
             distfun = function(c) as.dist(1 - c), axis=F, breaks=breaks, margins=c(0,0), 
             key = F, Rowv=NA, lower.left.only=FALSE, heatmapOnly=T, add.expr=abline(v=c(1.5,3.5,18.5,28.5))) 
