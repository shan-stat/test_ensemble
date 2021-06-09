# Required packages #
library(parallel); library(MASS); library(e1071); library(stats); library(boot); library(reticulate)
library(kyotil); library(aucm); library(cvAUC); library(dplyr); library(WeightedROC); library(glmnet); library(survey)
library(ranger); library(ebmc); library(randomForest); library(caret); library(SuperLearner)



# Simulation dataset in Section 2 #
sim.1=function(n1, n0, seed, sep=0) {
  set.seed(seed)
  Sigma <- matrix(c(1 , 0.7 , 0.5 , 0.3, 0.7 , 1 , 0.7 , 0.3, 0.5 , 0.7 , 1 , 0.3, 0.3 , 0.3 , 0.3 , 1 ),4,4)
  x=mvrnorm(n = n0+n1, rep(0, ncol(Sigma)), Sigma)
  colnames(x)="X"%.%1:ncol(Sigma)
  list(case=x[1:n1, ]+sep, control=x[n1+1:n0,])
}

# Simulated RV144 dataset in Section 3 #
sim.rv144=function(n1, n0, seed, alpha, betas, beta.z.1=0, beta.z.2=0, n=NULL, full=FALSE) {
  
  set.seed(seed)
  if (is.null(n)) {
    n=n1+n0
    bstrat=rep(1,n)
  } else {
    # two phase sampling stratum
    bstrat = rbern(n,prob=c(.15,.2,.2,.2,.25),generalized=T)
  }
  
  shared=rnorm(n,0,sd=sqrt(0.1)) # shared immune response, but not associated with outcome
  
  # group 1 
  p=1
  x.1=rnorm(n,0,sd=sqrt(0.6))
  grp.1=matrix(rnorm(p*n,0,sd=sqrt(0.1)), ncol=p)
  grp.1=shared+x.1+grp.1
  
  # group 2 
  p=2
  x.2=rnorm(n,0,sd=sqrt(0.6))
  grp.2=matrix(rnorm(p*n,0,sd=sqrt(0.3)), ncol=p)
  grp.2=shared+x.2+grp.2
  
  # group 3
  p=15
  # first simulate a latent variable, which is associated with outcome
  x.3=rnorm(n,0,sd=sqrt(0.6))
  # now add noise
  grp.3=matrix(rnorm(p*n,0,sd=sqrt(0.3)), ncol=p)
  # make some variables more noisy than others
  grp.3=t(t(grp.3) * c(rep(1,5), rep(10,5), rep(30,5)))
  grp.3=shared+x.3+grp.3
  grp.3=scale(grp.3)#otherwise, they will be too big after exp 
  # exponentiate group 3
  grp.3=exp((grp.3+10)*2)
  grp.3=scale(grp.3, center=FALSE) # center is FALSE so that we can take log if necessary
  
  # group 4
  p=10
  # first simulate a latent variable, which is associated with outcome
  x.4=rnorm(n,0,sd=sqrt(0.6))
  # now add noise
  grp.4=matrix(rnorm(p*n,0,sd=sqrt(0.3)), ncol=p)
  # make some variables more noisy than others
  grp.4=t(t(grp.4) * c(rep(1,6), rep(10,4)))
  grp.4=shared+x.4+grp.4
  grp.4=scale(grp.4)
  
  p=36; grp.5=rnorm(n,0,sd=sqrt(0.5)) + matrix(rnorm(p*n,0,sd=sqrt(0.4)), ncol=p) + shared
  p=16; grp.6=rnorm(n,0,sd=sqrt(0.5)) + matrix(rnorm(p*n,0,sd=sqrt(0.4)), ncol=p) + shared    
  p= 6; grp.7=rnorm(n,0,sd=sqrt(0.5)) + matrix(rnorm(p*n,0,sd=sqrt(0.4)), ncol=p) + shared
  p= 6; grp.8=rnorm(n,0,sd=sqrt(0.4)) + matrix(rnorm(p*n,0,sd=sqrt(0.5)), ncol=p) + shared
  p= 6; grp.9=rnorm(n,0,sd=sqrt(0.4)) + matrix(rnorm(p*n,0,sd=sqrt(0.5)), ncol=p) + shared
  p= 6;grp.10=rnorm(n,0,sd=sqrt(0.3)) + matrix(rnorm(p*n,0,sd=sqrt(0.6)), ncol=p) + shared
  p= 6;grp.11=rnorm(n,0,sd=sqrt(0.3)) + matrix(rnorm(p*n,0,sd=sqrt(0.6)), ncol=p) + shared
  p= 6;grp.12=rnorm(n,0,sd=sqrt(0.3)) + matrix(rnorm(p*n,0,sd=sqrt(0.6)), ncol=p) + shared
  p= 6;grp.13=rnorm(n,0,sd=sqrt(0.2)) + matrix(rnorm(p*n,0,sd=sqrt(0.7)), ncol=p) + shared
  p= 6;grp.14=rnorm(n,0,sd=sqrt(0.2)) + matrix(rnorm(p*n,0,sd=sqrt(0.7)), ncol=p) + shared
  p= 6;grp.15=rnorm(n,0,sd=sqrt(0.2)) + matrix(rnorm(p*n,0,sd=sqrt(0.7)), ncol=p) + shared
  p= 6;grp.16=rnorm(n,0,sd=sqrt(0.2)) + matrix(rnorm(p*n,0,sd=sqrt(0.7)), ncol=p) + shared
  p=10;grp.17=                          matrix(rnorm(p*n,0,sd=sqrt(0.9)), ncol=p) + shared
  
  if(beta.z.1==0 & beta.z.2==0) {
    # do this instead of generate random samples so that previous results without z1 z2 can be reproduced
    z.1=rep(0,n)
    z.2=rep(0,n)        
  } else {
    z.1=rbern(n,0.5)-0.5
    z.2=rnorm(n)
  }
  y=rbern(n, expit(alpha+beta.z.1*z.1+beta.z.2*z.2+betas[1]*x.1+betas[2]*x.2+betas[3]*x.3+betas[4]*x.4+betas[5]*x.2*x.4)); #mean(y)
  
  dat=data.frame(y, z.1, z.2, bstrat, grp.1, grp.2, grp.3, grp.4, grp.5, grp.6, grp.7, grp.8, grp.9, grp.10, grp.11, grp.12, grp.13, grp.14, grp.15, grp.16, grp.17)
  names(dat)=c("y","z1","z2", "bstrat", "x"%.%1:(ncol(dat)-4))
  
  # two-phase sampling 
  if (all(bstrat==1)) {
    dat$ph2=1
    dat$wt=1
  } else {
    tab=with(dat, table(bstrat, y))
    # sample all cases
    dat$ph2=ifelse(dat$y==1, 1, 0)
    # inverse sampling prob as wt
    dat$wt=1
    # fpc: the number of observations in each stratum of the population
    dat$fpc=1
    # sample 1:5 case:control ratio
    for (k in 1:max(dat$bstrat)) {
      dat$ph2[sample(which(dat$bstrat==k & dat$y==0), 5*tab[k,2])]=1 # stratified sampling without replacement
      dat$wt [       which(dat$bstrat==k & dat$y==0)]=tab[k,1]/(5*tab[k,2])
      dat$fpc[       which(dat$bstrat==k)]=table(bstrat)[k]
    }
  }
  
  # remove empty strata
  if(any(table(bstrat, y)[,'1']==0)){
    ept <- unname(which(table(bstrat, y)[,'1']==0))
    for(i in 1:length(ept)){
      dat <- subset(dat, bstrat!=ept[i])
    }
  }
  
  # phase two dataset vs full dataset
  if(full==FALSE){
    dat[dat$ph2==1,]
  } else if(full==TRUE){
    dat
  }
}

# k-fold CV #
kfold.split=function(k, n1, n0){
  training.subsets=list()
  test.subsets=list()
  tmp1=sample(1:n1)
  tmp0=sample(1:n0)
  splits=list()
  for (ki in 1:k) {
    splits[[ki]]=list(training=list(case=tmp1[(1:n1)%%k!=ki-1], control=tmp0[(1:n0)%%k!=ki-1]),
                      test=list(case=tmp1[(1:n1)%%k==ki-1], control=tmp0[(1:n0)%%k==ki-1]))
  }        
  splits
}

# Random 4:1 #
ran.kfold.split=function(k, n1, n0, replicates){
  training.subsets=list()
  test.subsets=list()
  splits=list()
  for (r in 1:replicates) {
    tmp1=sample(1:n1)
    tmp0=sample(1:n0)
    splits[[r]]= list(training=list(case=tmp1[(1:n1)%%k!=1], control=tmp0[(1:n0)%%k!=1]),
                      test=list(case=tmp1[(1:n1)%%k==1], control=tmp0[(1:n0)%%k==1]))
  }        
  splits
}

# Leave pair out CV #
lpo.split=function(n1, n0){
  training.subsets=list()
  test.subsets=list()
  splits=list()
  ind=0
  for (i in 1:n1) {
    for (j in 1:n0) {
      ind=ind+1
      splits[[ind]]=list(training=list(case=setdiff(1:n1,i), control=setdiff(1:n0,j)),
                         test=    list(case=i,               control=j))
    }        
  }        
  splits
}

# CV schemes #
get.splits=function(dat, cv.scheme=c("LPO","5fold","50xrandom4:1"), seed) {
  # save rng state before set.seed in order to restore before exiting this function
  save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
  if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }      
  set.seed(seed)    
  
  n0=nrow(dat$control)
  n1=nrow(dat$case)
  
  cv.scheme<-match.arg(cv.scheme)  
  if (cv.scheme=="LPO") {
    splits=lpo.split(n1, n0)
  } else if (cv.scheme=="5fold") {
    splits=kfold.split(5, n1, n0)
  } else if (cv.scheme=="50xrandom4:1") {
    splits=ran.kfold.split(5, n1, n0, 50)
  }
  # restore rng state 
  assign(".Random.seed", save.seed, .GlobalEnv)     
  splits
}

# Get cvAUC using multiple logistic regression #
get.cv.auc=function(dat, cv.scheme, seed=1){
  # k-fold CV #
  splits <- get.splits(dat, cv.scheme, seed)
  
  cv.aucs <-  mclapply( splits, function(split){
    dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
    dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
    set.seed(123)
    fit.mlr <- glm( factor(Y)~., dat.train, family=binomial ) 
    pred.mlr <- predict( fit.mlr, newdata=dat.test )
    fast.auc( pred.mlr, dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
  }, mc.cores = 4 )
  cv.aucs <- unlist( cv.aucs )
  mean(cv.aucs)
}

# Lasso screening #
screen_lasso <- function(Y, X, family, obsWeights=rep(1, nrow(X))) {
  set.seed(123)
  res.ls <- cv.glmnet(x=as.matrix(X), y=as.matrix(Y), weights=obsWeights, family=family, type.measure='auc', nfolds=5, alpha=1) # Lasso penalty
  vars.ls <- (coef(res.ls, s=res.ls$lambda.min) != 0)[-1]
  vars <- vars.ls; names(vars) <- colnames(X)
  return(vars)
}

# Stacking minimum P value and random forest #
get.st.auc = function(dat.train, mp.index, var.index, method, obsWeights=rep(1,nrow(dat.train)), seed=1){
  
  # 5-fold CV (fitting learners and obtaining out-of-sample prediction scores) #
  my_control <- trainControl(
    method="cv",
    number=5,
    savePredictions="final",
    classProbs=TRUE,
    summaryFunction=twoClassSummary
  )
  
  # Method of combining predictions from different learners #
  method <- get(method, mode = 'function')()
  if(!is.null(method$require)) {
    sapply(method$require, function(x) require(force(x), character.only = TRUE))
  }
  
  # Fitting candidate learners #
  # Target variable must be named as 'Y' #
  dat.train <- cbind(Y=dat.train$y, dat.train[,startsWith(names(dat.train),"z")|startsWith(names(dat.train),"x")])
  # Y=1,0 must be converted to factor that has two levels of 'case' (for Y=1) and 'control' (for Y=0) #
  dat.train$Y <- as.factor(dat.train$Y) ; levels(dat.train$Y) <- c('control','case')
  
  # Setting three RF hyper-parameters to their default values (without the following code, hyper-parameter tuning is automatically done) #
  rf.grid <- expand.grid(mtry = floor(sqrt(sum(var.index$rf)-1)), splitrule = 'gini', min.node.size = 1)
  set.seed( 123 )
  model_list <- caretList(
    list(Y~., data=dat.train[,var.index$rf], tuneGrid = rf.grid, weights = obsWeights), # RF with full features
    trControl=my_control,
    methodList=c('ranger')
  )
  
  # Meta-learner #
  set.seed( 123 )
  pred.cv <- makePredObsMatrix(model_list)
  pred.mat <- cbind(pred.cv$preds, "MP"=dat.train[,paste0("x",mp.index)])
  res.st.fit <- method$computeCoef(Z = pred.mat, Y =(as.numeric(pred.cv$obs)-1), obsWeights=obsWeights, libraryNames=names(var.index), verbose=FALSE, control=list(trimLogit=0.001))
  pred.st <- method$computePred(predY = pred.mat, coef = res.st.fit$coef, control=list(trimLogit=0.001))
  pred.st
}
