# Packages #
library( parallel ); library( foreach ); library( doParallel )
library( MASS ); library( e1071 ); library( stats )
library( kyotil ); library( aucm ); library( cvAUC )
library( glmnet ); library( lmtest )
library( rpart ); library( randomForest ); library( fastAdaboost )
library( ebmc )



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
        # sample 1:5 case:control ratio
        for (k in 1:max(dat$bstrat)) {
            dat$ph2[sample(which(dat$bstrat==k & dat$y==0), 5*tab[k,2])]=1
            dat$wt [       which(dat$bstrat==k & dat$y==0)]=tab[k,1]/(5*tab[k,2])
        }
    }
    
    dat[dat$ph2==1,]
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

# Getting CV-AUC (Likelihood-based function, machine learning) #
get.cv.auc <- function(dat, cv.scheme, method = c('MLR','RR','LR','DT','BG','RF','AB',
                                                  'UBG','OBG','SMOTEBG','URF','ORF','RUSB','SMOTEB'), numCores, seed=1){
  # Split data for k-fold CV #
  splits <- get.splits(dat, cv.scheme, seed)
  
  # MLR : Multiple Logistic Regression #
  if( method == 'MLR' ){
    cv.aucs <-  mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      set.seed(123)
      fit.mlr <- glm( factor(Y)~., dat.train, family=binomial ) 
      pred.mlr <- predict( fit.mlr, newdata=dat.test )
      fast.auc( pred.mlr, dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
    }, mc.cores = numCores )
    cv.aucs <- unlist( cv.aucs )
  }
  
  # RR : Ridge Regression #
  if( method == 'RR' ){
    cv.aucs <-  mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      dat.train.X <- as.matrix(dat.train[, colnames(dat.train) != 'Y']) ; dat.train.y <- as.matrix(dat.train[,'Y'])
      dat.test.X <- as.matrix(dat.test[,colnames(dat.test) != 'Y']) ; dat.test.y <- as.matrix(dat.test[,'Y'])
      set.seed(123)
      res.rr <- cv.glmnet( x = dat.train.X, y =  dat.train.y, family = "binomial" , type.measure = 'auc', nfolds = 5, alpha = 0 ) # Ridge penalty
      pred.rr <- as.numeric(drop( predict( res.rr , newx = dat.test.X , s = res.rr$lambda.min ) ))
      fast.auc( pred.rr, dat.test.y, reverse.sign.if.nece = FALSE, quiet = TRUE )
    }, mc.cores = numCores )
    cv.aucs <- unlist( cv.aucs )
  }
  
  # LR : Lasso Regression #
  if( method == 'LR' ){
    cv.aucs <-  mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      dat.train.X <- as.matrix(dat.train[, colnames(dat.train) != 'Y']) ; dat.train.y <- as.matrix(dat.train[,'Y'])
      dat.test.X <- as.matrix(dat.test[,colnames(dat.test) != 'Y']) ; dat.test.y <- as.matrix(dat.test[,'Y'])
      set.seed(123)
      res.lr <- cv.glmnet( x = dat.train.X, y =  dat.train.y, family = "binomial" , type.measure = 'auc', nfolds = 5, alpha = 1 ) # Lasso penalty
      pred.lr <- as.numeric(drop( predict( res.lr , newx = dat.test.X , s = res.lr$lambda.min ) ))
      fast.auc( pred.lr, dat.test.y, reverse.sign.if.nece = FALSE, quiet = TRUE )
    }, mc.cores = numCores )
    cv.aucs <- unlist( cv.aucs )
  }
  
  # DT : Decision Tree #
  if( method == 'DT' ){
    cv.aucs <- mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      set.seed(123)
      fit.dt <- rpart( factor(Y)~., dat.train )
      pred.dt <- predict( fit.dt, newdata=dat.test, type="prob" )
      fast.auc( pred.dt[,2], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
    }, mc.cores = numCores )
    cv.aucs <- unlist( cv.aucs )
  }
  
  # BG : Bagging #
  if( method == 'BG'){
    cv.aucs <- mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      set.seed(123)
      fit.bg <- randomForest( factor(Y)~., dat.train, mtry=( ncol(dat.train)-1 ) )
      pred.bg <- predict( fit.bg, newdata=dat.test, type="prob" )
      fast.auc( pred.bg[,2], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
    }, mc.cores = numCores )
    cv.aucs <- unlist( cv.aucs )
  }
  
  # RF : Random Forest #
  if( method == 'RF' ){
    cv.aucs <- mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      set.seed(123)
      fit.rf <- randomForest( factor(Y)~., dat.train )
      pred.rf <- predict( fit.rf, newdata=dat.test, type="prob" )
      fast.auc( pred.rf[,2], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
    }, mc.cores = numCores )
    cv.aucs <- unlist( cv.aucs )
  }
  
  # AB : AdaBoost.M2 #
  if( method == 'AB' ){
    cv.aucs <- mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      dat.train$Y <- factor( dat.train$Y )
      set.seed(123)
      fit.ab <- adam2( Y ~., data = dat.train, size = 100, alg = 'cart' )
      pred.ab <- predict( fit.ab, newdata=dat.test )
      fast.auc( pred.ab, dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
    }, mc.cores = numCores )
    cv.aucs <- unlist( cv.aucs )
  }
  
  # UBG : UnderBagging #
  if( method == 'UBG' ){
    cv.aucs <- mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      set.seed(123)
      fit.ubg <- randomForest( factor(Y)~., dat.train, mtry = ncol(dat.train)-1, strata = dat.train$Y, 
                               sampsize = rep(sum(dat.train$Y==1), 2) )
      pred.ubg <- predict( fit.ubg, newdata=dat.test, type="prob" )
      fast.auc( pred.ubg[,2], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
    }, mc.cores = numCores )
    cv.aucs <- unlist( cv.aucs )
  }
  
  # OBG : OverBagging #
  if( method == 'OBG' ){
    cv.aucs <- mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      set.seed(123)
      idx1 <- sample( 1:sum(dat.train$Y==1), sum(dat.train$Y==0), replace = TRUE ) # Oversampling minority class
      dat.temp <- list( case = subset(dat.train, Y==1, select=-Y), control = subset(dat.train, Y==0, select=-Y) )
      dat.train.mod <- rbind( data.frame(Y=1, dat.temp$case[idx1,,drop=F]),   data.frame(Y=0, dat.temp$control) )
      fit.obg <- randomForest( factor(Y)~., dat.train.mod, mtry = ncol(dat.train.mod)-1, strata = dat.train.mod$Y, 
                               sampsize = rep(sum(dat.train.mod$Y==1), 2) )
      pred.obg <- predict( fit.obg, newdata=dat.test, type="prob" )
      fast.auc( pred.obg[,2], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
    }, mc.cores = numCores )
    cv.aucs <- unlist( cv.aucs )
  }
  
  # SMOTEBG : SMOTEBagging #
  if( method == 'SMOTEBG' ){
    cv.aucs <- mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      dat.train$Y <- as.factor( dat.train$Y )
      set.seed(123)
      fit.smotebg <- sbag( Y ~ ., data = dat.train, size = 100, alg = 'cart' )
      pred.smotebg <- predict( fit.smotebg, newdata = dat.test )
      fast.auc( pred.smotebg, dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
    }, mc.cores = numCores )
    cv.aucs <- unlist( cv.aucs )
  }
  
  # URF : UnderRandomForest #
  if( method == 'URF' ){
    cv.aucs <- mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      set.seed(123)
      fit.urf <- randomForest( factor(Y)~., dat.train, strata = dat.train$Y, sampsize = rep(sum(dat.train$Y==1), 2) )
      pred.urf <- predict( fit.urf, newdata=dat.test, type="prob" )
      fast.auc( pred.urf[,2], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
    }, mc.cores = numCores )
    cv.aucs <- unlist( cv.aucs )
  }
  
  # ORF : OverRandomForest #
  if( method == 'ORF' ){
    cv.aucs <- mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      
      # Oversampling minority class
      set.seed(123)
      idx1 <- sample( 1:sum(dat.train$Y==1), sum(dat.train$Y==0), replace = TRUE ) 
      dat.temp <- list( case = subset(dat.train, Y==1, select=-Y), control = subset(dat.train, Y==0, select=-Y) )
      dat.train.mod <- rbind( data.frame(Y=1, dat.temp$case[idx1,,drop=F]),   data.frame(Y=0, dat.temp$control) )
      
      set.seed(123)
      fit.orf <- randomForest( factor(Y)~., dat.train.mod, strata = dat.train.mod$Y, sampsize = rep(sum(dat.train.mod$Y==1), 2) )
      pred.orf <- predict( fit.orf, newdata=dat.test, type="prob" )
      fast.auc( pred.orf[,2], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
    }, mc.cores = numCores )
    cv.aucs <- unlist( cv.aucs )
  }
  
  # RUSB : Random Under Sampling Boosting #
  if( method == 'RUSB' ){
    cv.aucs <- mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      dat.train$Y <- as.factor( dat.train$Y )
      set.seed(123)
      fit.rusb <- rus( Y ~ ., data = dat.train, size = 100, alg = 'cart', ir = 1 )
      pred.rusb <- predict( fit.rusb, newdata = dat.test )
      fast.auc( pred.rusb, dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
    }, mc.cores = numCores )
    cv.aucs <- unlist( cv.aucs )
  }
  
  # SMOTEB : SMOTEBoost #
  if( method == 'SMOTEB' ){
    cv.aucs <- mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      dat.train$Y <- as.factor( dat.train$Y )
      set.seed(123)
      fit.smoteb <- sbo( Y ~ ., data = dat.train, size = 100, alg = 'cart', over = 500 )
      pred.smoteb <- predict( fit.smoteb, newdata = dat.test )
      fast.auc( pred.smoteb, dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
    }, mc.cores = numCores )
    cv.aucs <- unlist( cv.aucs )
  }
  
  mean(cv.aucs)
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

# Wald test / LRT #
get.pvalue.test <- function( dat ){
  set.seed( 123 )
  null.model <- glm( factor(y) ~ 1, dat, family=binomial ) # without covariates
  
  # Generating result vector #
  pvalues.wald <- c()
  pvalues.lrt <- c()
  for( i in 1:(ncol(dat)-1) ){
    dat.temp <- dat[, c('y',paste('x',i,sep=''))] # without covariates
    set.seed(123)
    fit.model <- glm( factor(y) ~ ., dat.temp, family=binomial ) 
    res.wald <- waldtest( null.model, fit.model )
    pvalues.wald[i] <- res.wald$'Pr(>F)'[2]
    res.lrt <- lrtest( null.model, fit.model )
    pvalues.lrt[i] <- res.lrt$'Pr(>Chisq)'[2]
  }
    
  c(min(pvalues.wald), min(pvalues.lrt))
}

# CV-AUC for Univariate Logistic Regression (for Section 3.1)#
get.cv.auc.ULR <- function( dat, cv.scheme = c('5fold','nocv'), numCores = numCores, seed=1 ){
  # CV-AUC #
  if( cv.scheme == '5fold' ){
    # Split data for k-fold CV #
    splits <- get.splits(dat, cv.scheme, seed = seed)
    # Generating result vector #
    cvaucs <- c()
    for( i in 1:ncol(dat$case) ){
      cv.aucs <-  mclapply( splits, function(split){
        dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case, paste('x',i,sep=''),drop=F]),   data.frame(Y=0, dat$control[split$training$control, paste('x',i,sep=''),drop=F]) )
        dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case, paste('x',i,sep=''),drop=F]),       data.frame(Y=0, dat$control[split$test$control, paste('x',i,sep=''),drop=F]) )
        set.seed(123)
        fit.lr <- glm( factor(Y) ~ ., dat.train, family=binomial ) 
        pred.lr <- predict( fit.lr, newdata=dat.test )
        fast.auc( pred.lr, dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
      }, mc.cores = numCores )
      cvaucs[i] <- mean( unlist( cv.aucs ) )
    }
  } 
  
  # AUC #
  if( cv.scheme == 'nocv' ){
    # Generating result vector #
    cvaucs <- c()
    for( i in 1:ncol(dat$case) ){
      dat.train <- rbind( data.frame(Y=1, dat$case[,c(paste('x',i,sep='')),drop=F]), data.frame(Y=0, dat$control[,c(paste('x',i,sep='')),drop=F]) )
      set.seed(123)
      fit.lr <- glm( factor(Y) ~ ., dat.train, family=binomial ) 
      pred.lr <- predict( fit.lr, newdata=dat.train )
      cvaucs[i] <- fast.auc( pred.lr, dat.train$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
    }
  }
    
  max(cvaucs)
}


lgb.normalizedgini = function(preds, dtrain){
  actual = getinfo(dtrain, "label")
  score  = NormalizedGini(preds,actual)
  return(list(name = "gini", value = score, higher_better = TRUE))
}                
