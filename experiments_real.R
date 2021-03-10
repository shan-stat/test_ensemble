### Import helper functions ###
setwd("set your local directory where helper_cvauc.R is located")
source("helper_cvauc.R")

# The following four files are required for fitting stacking #
setwd("set your local directory where the following four files are located")
source('helper_functions.R')
source('caretList.R')
source('caretEnsemble.R')
source('caretStack.R')



### 1. Real RV144 dataset ###
#setwd("/Users/shan/Desktop/Paper/YFong/1.ML")
#devtools::install_local("RV144cc_2019.08-14.tar")
data("allcc", package="RV144cc")
tmp=names(allcc)
tmp=tmp[!endsWith(tmp, "_cat")]
tmp=tmp[!startsWith(tmp,"c")]

# For using all markers #
imm.var.all = tmp[2:195] # using all markers
dat.all <- select(allcc, c("y","gender","behavioral_risk",imm.var.all))
names(dat.all) <- c("y","z1","z2","x"%.%1:194) # for simplicity, convert names to z and x.
dat <- data.frame(dat.all, ipw=1/allcc$sampling.p)

# For using all markers without IgAC1 #
imm.var = tmp[2:195]
imm.var.igac1 = imm.var[-49] # using all markers without IgAC1
dat.igac1 <- select(allcc, c("y","gender","behavioral_risk",imm.var.igac1))
names(dat.igac1) <- c("y","z1","z2","x"%.%1:193) # for simplicity, convert names to z and x.
dat <- data.frame(dat.igac1, ipw=1/allcc$sampling.p)

# Mean imputation for missing values #
for(i in 1:ncol(dat)){
  if(sum(is.na(dat[,i])) != 0){
    set.seed(123)
    dat[which(is.na(dat[,i])),i] <- mean(dat[,i], na.rm=TRUE) # mean
  }
}



### 2. Setting process arguments ###
Args <- commandArgs(trailingOnly=TRUE) 
if (length(Args)==0) {
  # batch.size and batch.number are for a batch script
  # sim.setting: rv144_1 for real RV144 experiments
  # fit.setting: 5fold for 5-fold cross validation
  # proj: you can choose one of testing methods
  # ipw: uw for unweighted; sw for semi-weighted; w for weighted
  # scr: TRUE for lasso screening; FALSE for no screening
  # stacking: TRUE for stacking (note that to use stacking, scr should be set to scr="TRUE")
  Args=c(batch.size="1",batch.number="1",setting="rv144_1",fit.setting="5fold",proj="perm_RF",scr="TRUE",stacking="FALSE")
}
myprint(Args)
i=0;
i=i+1; batch.size=as.numeric(Args[i])
i=i+1; batch=as.numeric(Args[i])
seeds=1:batch.size+batch.size*(batch-1); names(seeds)=seeds
myprint(batch, batch.size)
i=i+1; setting=Args[i]; tmp=strsplit(setting,"_")[[1]]
model=tmp[1]
i=i+1; fit.setting=Args[i]
cv.scheme=fit.setting
i=i+1; proj=Args[i]
i=i+1; scr=Args[i]
i=i+1; stacking=Args[i]

nperm=1e3 # permutation replicates
verbose=ifelse(unix(),0,2)



### 3. Experiments ###
res=sapply(seeds, simplify="array", function (seed) {
  
  myprint(seed)
  
  # Real RV144 dataset #
  N=nrow(dat)
  n1=sum(dat$y)
  n0=N-n1
  p=sum(startsWith(names(dat),"x"))
  myprint(n1,n0,p)
  f="as.factor(Y)~z1+z2+" %.% concatList(names(dat)[startsWith(names(dat),"x")],"+")    
  
  # 2-1. Fitting testing methods #
  if (startsWith(proj,"BH")) {
    # BH: Benjamini-Hochberg #
    pvals=sapply (1:p, function(i) {
      set.seed(123)
      fit=glm(as.formula("y~z1+z2++x"%.%i), dat, family=binomial) # unweighted
      last(summary(fit)$coef)
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
          fit=glm(as.formula("y~z1+z2+x"%.%i), dat.train, family=binomial) # unweighted
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
            dat.train <- select(dat.train,c("Y","y","z1","z2",names(which(screened.var)),"ipw"))
            dat.test <- select(dat.test,c("Y","y","z1","z2",names(which(screened.var)),"ipw"))
            # convert factor to dummy variable because cv.glmnet cannot deal with factor directly #
            z1.train.dummy <- model.matrix(~z1-1, dat.train); dat.train$z1 <- z1.train.dummy
            z2.train.dummy <- model.matrix(~z2-1, dat.train); dat.train$z2 <- z2.train.dummy 
            z1.test.dummy <- model.matrix(~z1-1, dat.test); dat.test$z1 <- z1.test.dummy
            z2.test.dummy <- model.matrix(~z2-1, dat.test); dat.test$z2 <- z2.test.dummy 
            dat.train.X <- as.matrix(select( dat.train, -c('Y','y','ipw') )); dat.train.y <- as.matrix(select( dat.train, 'Y' ))
            dat.test.X <- as.matrix(select( dat.test, -c('Y','y','ipw') )); dat.test.y <- as.matrix(select( dat.test, 'Y' ))
            set.seed(123)
            fit.rr <- cv.glmnet(x=dat.train.X, y=dat.train.y, family="binomial", type.measure='auc', nfolds=5, alpha=0) # cv.glmnet cannot directly deal with factor (dummy variables are needed)
            pred.rr <- as.numeric(drop( predict( fit.rr, newx=dat.test.X, s=fit.rr$lambda.min ) ))
            fast.auc( pred.rr, dat.test.y, reverse.sign.if.nece = FALSE, quiet = TRUE )
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
            fit.bg <- randomForest(as.formula(f), data=dat.train, mtry=(sum(screened.var)+2))
            pred.bg <- predict(fit.bg, newdata=dat.test, type='prob')
            fast.auc( pred.bg[,2], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
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
          if(scr=="TRUE" | stacking=="TRUE"){
            vars <- names(dat.train)[startsWith(names(dat.train),"x")]
            screened.var <- screen_lasso(Y=dat.train$Y, X=dat.train[,vars], family="binomial") # lasso
          }
          
          # Fitting #
          if(scr=="TRUE" & sum(screened.var)==0){
            0.5
          } else {
            if(stacking=="TRUE"){
              # For only stacking #
              all.var <- rep(TRUE,sum(startsWith(names(dat.train),"x")))
              var.index.set <- list(rf=c(TRUE,TRUE,TRUE,all.var), rf.scr=c(TRUE,TRUE,TRUE,screened.var)) # the first three TRUE are for Y, z1, and z2
              set.seed(123)
              pred.st <- get.rf.st.cvauc(dat.train=dat.train, dat.test=dat.test, var.index=var.index.set, method='method.NNloglik')
              fast.auc( pred.st, dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
            } else {
              # For no screening and screening #
              f = "as.factor(Y)~z1+z2+" %.% concatList(names(which(screened.var)),"+")
              set.seed(123)
              fit.rf <- randomForest(as.formula(f), data=dat.train)
              pred.rf <- predict(fit.rf, newdata=dat.test, type='prob')
              fast.auc( pred.rf[,2], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
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
            fit.ab <- adam2(Y~., data=select(dat.train, c("Y","z1","z2",names(which(screened.var)))), size=100, alg='cart')
            pred.ab <- predict( fit.ab, newdata=dat.test )
            fast.auc( pred.ab, dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
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
            fit.ubg <- ranger( as.formula(f), data=dat.train, case.weights=weights.samp, probability=TRUE, min.node.size=1, 
                               mtry=sum(screened.var)+2, sample.fraction = (class.dist['1']*2)/sum(class.dist))
            pred.ubg <- predict( fit.ubg, data=dat.test )
            fast.auc( pred.ubg$predictions[,'1'], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
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
            fit.urf <- ranger( as.formula(f), data=dat.train, case.weights=weights.samp, probability=TRUE, min.node.size=1, sample.fraction = (class.dist['1']*2)/sum(class.dist) )
            pred.urf <- predict( fit.urf, data=dat.test )
            fast.auc( pred.urf$predictions[,'1'], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
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
            fit.rusb <- rus(Y~., data=select(dat.train, c("Y","z1","z2",names(which(screened.var)))), size=100, alg='cart')
            pred.rusb <- predict( fit.rusb, newdata = dat.test )
            fast.auc( pred.rusb, dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
          }
        })
        mean(cv.aucs)
        
      } else if(startsWith(proj,"perm_overBG")) {
        # overBG: bagging with over-sampling #
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
          
          # over-sampling minority class #
          set.seed(123)
          idx1 <- sample( 1:sum(dat.train$Y==1), sum(dat.train$Y==0), replace = TRUE )
          dat.train.over <- rbind(dat.train[idx1,,drop=F], subset(dat.train, Y==0))
          
          # Fitting #
          if(scr=="TRUE" & sum(screened.var)==0){
            0.5
          } else {
            f = "as.factor(Y)~z1+z2+" %.% concatList(names(which(screened.var)),"+")
            set.seed(123)
            fit.obg <- ranger( as.formula(f), data=dat.train.over, probability=TRUE, min.node.size=1, mtry=sum(screened.var)+2)
            pred.obg <- predict( fit.obg, data=dat.test )
            fast.auc( pred.obg$predictions[,'1'], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
          }
        })
        mean(cv.aucs)
        
      } else if(startsWith(proj,"perm_overRF")) {
        # underRF: random forest with over-sampling #
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
          
          # over-sampling minority class #
          set.seed(123)
          idx1 <- sample( 1:sum(dat.train$Y==1), sum(dat.train$Y==0), replace = TRUE )
          dat.train.over <- rbind(dat.train[idx1,,drop=F], subset(dat.train, Y==0))
          
          # Fitting #
          if(scr=="TRUE" & sum(screened.var)==0){
            0.5
          } else {
            f = "as.factor(Y)~z1+z2+" %.% concatList(names(which(screened.var)),"+")
            set.seed(123)
            fit.orf <- ranger( as.formula(f), data=dat.train.over, probability=TRUE, min.node.size=1 )
            pred.orf <- predict( fit.orf, data=dat.test )
            fast.auc( pred.orf$predictions[,'1'], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
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
            fit.smoteb <- sbo( Y~., data=select(dat.train, c("Y","z1","z2",names(which(screened.var)))), size=100, alg='cart', over=500 )
            pred.smoteb <- predict( fit.smoteb, newdata = dat.test )
            fast.auc( pred.smoteb, dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
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
      ind.1=sample(1:N, n1)    # cases in the permuted dataset
      ind.0=setdiff(1:N,ind.1) # controls in the permuted dataset
      dat.b=list()
      dat.b$case=   cbind(dat.case.control[1:n1,   c("z1","z2","ipw")], subset(dat.case.control[ind.1,], select=-c(z1,z2,ipw)) )
      dat.b$control=cbind(dat.case.control[1:n0+n1,c("z1","z2","ipw")], subset(dat.case.control[ind.0,], select=-c(z1,z2,ipw)) )
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
