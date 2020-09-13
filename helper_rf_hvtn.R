# Packages #
library(parallel); library(foreach); library(doParallel)
library(aucm); library(cvAUC) ; library(MASS)
library(randomForest)
library(HVTN505) ; library(kyotil) ; library(dplyr)


# K-fold CV #
kfold.split=function(k, n1, n0, seed){
  training.subsets=list()
  test.subsets=list()
  set.seed(seed)
  tmp1=sample(1:n1)
  tmp0=sample(1:n0)
  splits=list()
  for (ki in 1:k) {
    splits[[ki]]=list(training=list(case=tmp1[(1:n1)%%k!=ki-1], control=tmp0[(1:n0)%%k!=ki-1]),
                      test=list(case=tmp1[(1:n1)%%k==ki-1], control=tmp0[(1:n0)%%k==ki-1]))
  }        
  splits
}

# Split function #
get.splits=function(dat, cv.scheme=c("5fold"), seed=1) {
  set.seed(seed)
  if (cv.scheme=='5fold') {
    n1=nrow(dat$case)
    n0=nrow(dat$control)
    splits=kfold.split(5, n1, n0, seed)
  }
  splits
}

# Variable names #
get_nms_group_all_antigens <- function(X, assays, assays_to_exclude = "") {
  ## set all vars to be false
  vars <- rep(FALSE, ncol(X))
  ## set vars with assay in name to be true
  ## may be more than one
  for (i in 1:length(assays)) {
    if (assays_to_exclude != "") {
      vars[grepl(assays[i], names(X)) & !grepl(assays_to_exclude, names(X))] <- TRUE
    } else {
      vars[grepl(assays[i], names(X))] <- TRUE
      if (assays[i] == "phago") {
        vars[grepl("ADCP1", names(X))] <- TRUE
      }
    }
    
  }
  return(vars)
}

# CV-AUC (k-fold CV) #
get.kfold.cvauc = function(dat, cv.scheme, method=c('RF'), seed=1){
  
  # Split data for k-fold CV #
  splits <- get.splits(dat, cv.scheme, seed)
  
  # Random Forest #
  if( method == 'RF' ){
    cv.aucs <-  mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      set.seed( 123 )
    
    

# GLM
seeds=1:10; names(seeds)=seeds
out=sapply (seeds, simplify="array", function(seed){        
      splits <- get.splits(dat, cv.scheme, seed)
      folds=1:5; names(folds)=folds
      sapply(folds, simplify="array", function(f) {
            split=splits[[f]]
            dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
            dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
    
            fit <- glm( Y~., data = dat.train,family="binomial")
            
            # train, to see if it overfits
            trai.rf <- predict( fit, newdata=dat.train, type = 'resp' )    
            pred.rf <- predict( fit, newdata=dat.test, type = 'resp' )
            c(train=fast.auc( trai.rf, dat.train$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
            ,predi=fast.auc( pred.rf, dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE ))
       })
})
glm.perf=apply(out, 1, mean)







    # 
    params=list(
        c(mtry=3, maxnodes=2),
        c(mtry=3, maxnodes=3),
        c(mtry=3, maxnodes=4),
        c(mtry=3, maxnodes=10)
    )
    
    tabs=list()
    for (i in 1:4) {
        seeds=1:10; names(seeds)=seeds
        out=sapply (seeds, simplify="array", function(seed){        
              splits <- get.splits(dat, cv.scheme, seed)
              folds=1:5; names(folds)=folds
              sapply(folds, simplify="array", function(f) {
                    split=splits[[f]]
                    dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
                    dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
                    ntrees=seq(10,500,by=50)
                    out=sapply (ntrees, function(ntree) {        
                        fit.rf <- randomForest( factor(Y)~., data = dat.train, mtry=params[[i]]["mtry"], nodesize=5, maxnodes=params[[i]]["maxnodes"], ntree=ntree)
                        
                        # train, to see if it overfits
                        trai.rf <- predict( fit.rf, newdata=dat.train, type = 'prob' )    
                        pred.rf <- predict( fit.rf, newdata=dat.test, type = 'prob' )
                        c(train=fast.auc( trai.rf[,2], dat.train$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
                        ,predi=fast.auc( pred.rf[,2], dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE ))
                    })              
               })
        })
        tab=apply(out, 1:2, mean)
        tabs[[i]]=tab
    }

    myfigure(mfrow=c(1,4))
    for (i in 1:4) {
        mymatplot(ntrees, t(tabs[[i]]), ylim=c(0, 1), xlab="ntree", main=paste0("mtry=",params[[i]]["mtry"], ", maxnodes=",params[[i]]["maxnodes"]))
        abline(h=.5,col="darkgray")
        abline(h=glm.perf, lty=c(1:2))
    }



      library(reprtree) # download from github and remove tree.rd and compile
      reprtree:::plot.getTree(fit.rf)
        
      myfigure(mfrow=c(3,1))
      reprtree:::plot.getTree(tr=reprtree:::as.tree(getTree(fit.rf, 1, labelVar=TRUE),fit.rf))
      reprtree:::plot.getTree(tr=reprtree:::as.tree(getTree(fit.rf, 2, labelVar=TRUE),fit.rf))
      reprtree:::plot.getTree(tr=reprtree:::as.tree(getTree(fit.rf, 3, labelVar=TRUE),fit.rf))



              


    }, mc.cores = 4 )
    cv.aucs <- unlist( cv.aucs )
  }
  cv.aucs
}


fit=glm(Y~., dat.train, family="binomial")
summary(fit)
      pred.rf <- predict( fit, newdata=dat.test, type = 'resp' )
      fast.auc( pred.rf, dat.test$Y, reverse.sign.if.nece = FALSE, quiet = TRUE )
