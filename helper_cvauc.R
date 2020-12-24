# Required packages #
library(parallel); library(foreach); library(doParallel)
library(MASS); library(e1071); library(stats)
library(kyotil); library(aucm); library(cvAUC)
library(glmnet); library(lmtest)
library(rpart); library(randomForest); library(fastAdaboost); library(ebmc)



# simulate a dataset to mimic the RV144 dataset with 150 variables 
# Signals: 
# Group 1. IgG V2 binding
# Group 2. IgA binding
# Group 3. T-cell cytokines. The signal is on the log scale, but the measured variable is on the linear scale
# Group 4. Avidity
# 13 more groups
# all variables share at least 0.1 similarity
sim.rv144=function(n1, n0, seed, alpha, betas, beta.z.1=0, beta.z.2=0, n=NULL) {
  
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
    p=7
    # first simulate a latent variable, which is associated with outcome
    x.1=rnorm(n,0,sd=sqrt(0.6))
    # now add noise
    grp.1=matrix(rnorm(p*n,0,sd=sqrt(0.3)), ncol=p)
    # make some variables more noisy than others
    grp.1=t(t(grp.1) * c(1, rep(100,6)))
    grp.1=shared+x.1+grp.1
    grp.1=scale(grp.1)
    
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
    
    p=30; grp.5=rnorm(n,0,sd=sqrt(0.5)) + matrix(rnorm(p*n,0,sd=sqrt(0.4)), ncol=p) + shared
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
    
    #alpha=-1.9 # adjusted so that the prevalence is 1/6 with the following beta's
    #beta.1=2/3; beta.2=2/3; beta.12=-2/3; beta.3=2/3
    
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
            dat$ph2[sample(which(dat$bstrat==k & dat$y==0), 5*tab[k,2])]=1
            dat$wt [       which(dat$bstrat==k & dat$y==0)]=tab[k,1]/(5*tab[k,2])
            dat$fpc[       which(dat$bstrat==k)]=table(bstrat)[k]
        }
    }
  
    # remove empty strata for svyglm
    if(any(table(bstrat, y)[,'1']==0)){
      ept <- unname(which(table(bstrat, y)[,'1']==0))
      for(i in 1:length(ept)){
        dat <- subset(dat, bstrat!=ept[i])
      }
    }
  
    # for stratified independent sampling without replacement
    dat[dat$ph2==1,]
    # for phase-two sampling scheme
    #dat
}

#seed=1; 
#alpha=-5; betas=rep(0,5) # 717 cases for n=1e5
##alpha=-6.6; betas=c(0.9,0.7,1.2,-1,-1.2); # 731 cases for n=1e5
#beta.z.1=0.5; beta.z.2=0
#dat=sim.rv144(n1=NULL, n0=NULL, seed, alpha, betas, beta.z.1, beta.z.2, n=1e5) 
#with(dat, table(bstrat, ph2, y))
#with(dat, table(wt, ph2, y))
#sum(dat$y)



# simulate a dataset to mimic the MTCT dataset
# 50 variables
# 1 signal
sim.mtct=function(n1, n0, seed, alpha, beta, beta.z.1=0, beta.z.2=0) {
  set.seed(seed)
  n=n1+n0
  # group 1:  5 variables, around 0.7 covariance
  # group 2: 10 variables, around 0.7 covariance
  # group 3: 10 variables, around 0.5 covariance
  # group 4: 10 variables, around 0.3 covariance
  # group 5: 15 variables, independent
  
  # group 1
  # first simulate a latent variable
  p=5
  x.1=rnorm(n,0,sd=sqrt(0.7)) # x.1 will be associated with outcome
  grp.1=x.1+matrix(rnorm(p*n,0,sd=sqrt(0.3)), ncol=p)
  
  p=10; grp.2=rnorm(n,0,sd=sqrt(0.7)) + matrix(rnorm(p*n,0,sd=sqrt(0.3)), ncol=p)
  p=10; grp.3=rnorm(n,0,sd=sqrt(0.5)) + matrix(rnorm(p*n,0,sd=sqrt(0.5)), ncol=p)
  p=10; grp.4=rnorm(n,0,sd=sqrt(0.3)) + matrix(rnorm(p*n,0,sd=sqrt(0.7)), ncol=p)
  p=15; grp.5=                          matrix(rnorm(p*n,0,sd=sqrt(1  )), ncol=p)
  
  #beta=0.5; beta.z.1=2; beta.z.2=2
  if(beta.z.1==0 & beta.z.2==0) {
    # do this instead of generate random samples so that previous results without z1 z2 can be reproduced
    z.1=rep(0,n)
    z.2=rep(0,n)        
  } else {
    z.1=rbern(n,0.5)-0.5
    z.2=rnorm(n)
  }
  #alpha=-0.7
  y=rbern(n, expit(alpha+beta.z.1*z.1+beta.z.2*z.2+beta*x.1)); mean(y)
  
  dat=data.frame(y, z.1, z.2, grp.1, grp.2, grp.3, grp.4, grp.5)
  names(dat)=c("y","z1","z2","x"%.%1:50)
  dat
}


sim.1=function(n1, n0, seed, sep=0) {
    set.seed(seed)
    Sigma <- matrix(c(1 , 0.7 , 0.5 , 0.3, 0.7 , 1 , 0.7 , 0.3, 0.5 , 0.7 , 1 , 0.3, 0.3 , 0.3 , 0.3 , 1 ),4,4)
    x=mvrnorm(n = n0+n1, rep(0, ncol(Sigma)), Sigma)
    colnames(x)="X"%.%1:ncol(Sigma)
    list(case=x[1:n1, ]+sep, control=x[n1+1:n0,])
}


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

# replicates is the number of splits needed
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
#ran.kfold.split(5,10,20,3)

# leave (one) pair out
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
#lopo.split(3,4)

# CV schemes #
get.splits=function(dat, cv.scheme=c("5fold","random5fold","LPO","nocv"), seed) {
    # save rng state before set.seed in order to restore before exiting this function
    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
    if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }      
    set.seed(seed)    
    
    n0=nrow(dat$control)
    n1=nrow(dat$case)
    
    cv.scheme<-match.arg(cv.scheme)  
    if (cv.scheme=="5fold") {
        splits=kfold.split(5, n1, n0)
    } else if (cv.scheme=="random5fold") {
        splits=ran.kfold.split(5, n1, n0, 50)
    } else if (cv.scheme=="LPO") {
        splits=lpo.split(n1, n0)
    } else if (cv.scheme=="nocv") {
        splits=list()
        splits[[1]]=list(training=list(case=(1:n1), control=(1:n0)),
                         test=list(case=(1:n1), control=(1:n0)))
    }
    # restore rng state 
    assign(".Random.seed", save.seed, .GlobalEnv)     
    splits
}

# LeDell's CI #
get.cv.auc.LeDell=function(dat, cv.scheme, seed) {
    splits=get.splits(dat, cv.scheme, seed)
    scores=lapply (splits, function (split){
        dat.train=rbind(data.frame(Y=1,dat$case[split$training$case,,drop=F]),   data.frame(Y=0,dat$control[split$training$control,,drop=F]))
        dat.test =rbind(data.frame(Y=1,dat$case[split$test$case,,drop=F]),       data.frame(Y=0,dat$control[split$test$control,,drop=F]))
        fit=glm(Y~., dat.train, family=binomial)    
        predict(fit, newdata=dat.test)
    })
    labels=lapply (splits, function (split){
        dat.test =rbind(data.frame(Y=1,dat$case[split$test$case,,drop=F]),       data.frame(Y=0,dat$control[split$test$control,,drop=F]))
        dat.test$Y
    })
    out <- ci.cvAUC(scores, labels)    
    c(est=out$cvAUC, lb=out$ci[1], ub=out$ci[2])
}

lgb.normalizedgini = function(preds, dtrain){
  actual = getinfo(dtrain, "label")
  score  = NormalizedGini(preds,actual)
  return(list(name = "gini", value = score, higher_better = TRUE))
}
