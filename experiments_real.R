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
#devtools::install_local("RV144cc_2019.08-14.tar")
data("allcc", package="RV144cc")
tmp=names(allcc)
tmp=tmp[!endsWith(tmp, "_cat")]; tmp=tmp[!startsWith(tmp,"c")] # without categorical variables

scenario <- c("full","-IgAC1")[1] # full or -IgAC1 
if(scenario=="full"){
  # Full scenario #
  imm.var.all = tmp[2:195] # using all markers
  dat.all <- select(allcc, c("y","gender","behavioral_risk",imm.var.all))
  names(dat.all) <- c("y","z1","z2","x"%.%1:194) # for simplicity, convert names to z (clinical covariates) and x (biomarkers).
  dat <- data.frame(dat.all, ipw=1/allcc$sampling.p)
} else if(scenario=="-IgAC1"){
  # -IgAC1 scenario #
  imm.var = tmp[2:195]
  imm.var.igac1 = imm.var[-49] # using all markers without IgAC1
  dat.igac1 <- select(allcc, c("y","gender","behavioral_risk",imm.var.igac1))
  names(dat.igac1) <- c("y","z1","z2","x"%.%1:193) # for simplicity, convert names to z (clinical covariates) and x (biomarkers).
  dat <- data.frame(dat.igac1, ipw=1/allcc$sampling.p)
}

# Mean imputation for missing values #
for(i in 1:ncol(dat)){
  if(sum(is.na(dat[,i])) != 0){
    set.seed(123)
    dat[which(is.na(dat[,i])),i] <- mean(dat[,i], na.rm=TRUE) # mean
  }
}



### 2. Setting arguments ###
Args <- commandArgs(trailingOnly=TRUE) 
if (length(Args)==0) {
  # batch.size and batch.number are for runnning a batch script
  # setting: "rv144_1" for real data experiment in Section 4
  # fit.setting: "5-fold" by default
  # proj: "perm_MP", "perm_RR", "perm_BG", "perm_RF", "perm_AB", "perm_ST"
  Args=c(batch.size="1",batch.number="1",setting="rv144_1",fit.setting="5fold",proj="perm_MP")
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
  if (startsWith(proj,"perm")) {   
    # perm: Permutation-based tests #
    do.est=function(dat.b){
      
      if(startsWith(proj,"perm_MP")) {
        # MP: the minimum p-value #
        dat.train=rbind(data.frame(y=1,dat.b$case), data.frame(y=0,dat.b$control))
        pvals=sapply (1:p, function(i) {
          set.seed(123)
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
        })
        mean(cv.aucs)
        
      } else if(startsWith(proj,"perm_BG")) {
        # BG: bagging #
        splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)
        cv.aucs <-  sapply( splits, function(split){
          dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
          dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
          set.seed(123)
          fit.bg <- randomForest(as.formula(f), data=dat.train, mtry=(length(names(dat)[startsWith(names(dat),"x")])+2))
          pred.bg <- predict(fit.bg, newdata=dat.test, type='prob')
          fast.auc( pred.bg[,2], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
        })
        mean(cv.aucs)
        
      } else if(startsWith(proj,"perm_RF")) {
        # RF: random forest #
        splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)
        cv.aucs <-  sapply( splits, function(split){
          dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
          dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
          set.seed(123)
          fit.rf <- randomForest(as.formula(f), data=dat.train)
          pred.rf <- predict(fit.rf, newdata=dat.test, type='prob')
          fast.auc( pred.rf[,2], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
        })
        mean(cv.aucs)
    
      } else if(startsWith(proj,"perm_AB")) {
        # AB: adaptive boosting #
        splits <- get.splits(dat.b, cv.scheme=cv.scheme, seed=1)
        cv.aucs <-  sapply( splits, function(split){
          dat.train <- rbind( data.frame(Y=1, dat.b$case[split$training$case,,drop=F]),   data.frame(Y=0, dat.b$control[split$training$control,,drop=F]) )
          dat.test <-  rbind( data.frame(Y=1, dat.b$case[split$test$case,,drop=F]),       data.frame(Y=0, dat.b$control[split$test$control,,drop=F]) )
          dat.train$Y <- factor(dat.train$Y)
          set.seed(123)
          fit.ab <- adam2(Y~., data=select(dat.train, c("Y","z1","z2",names(dat)[startsWith(names(dat),"x")])), size=100, alg='cart')
          pred.ab <- predict( fit.ab, newdata=dat.test )
          fast.auc( pred.ab, dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
        })
        mean(cv.aucs)
        
      } else if(startsWith(proj,"perm_ST")) {
        # ST: stacking minimum p-value and random forest #
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
