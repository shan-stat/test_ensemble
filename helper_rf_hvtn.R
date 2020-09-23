# Packages #
library(parallel); library(foreach); library(doParallel)
library(aucm); library(cvAUC); library(MASS)
library(HVTN505); library(kyotil); library(dplyr); library(vimp)
library(glmnet); library(tuneRanger); library(caret)
library(nnls); library(quadprog); library(nloptr)



# Logistic regression univariate p-value < level #
rank_univariate_logistic_pval_plus_exposure <- function(Y, X, family, obsWeights, ...) {
  # logistic regression of outcome on each variable
  listp <- apply(X, 2, function(x, Y, family) {
    summ <- coef(summary(glm(Y ~ x + X$age + X$BMI + X$bhvrisk, family = family, weights = obsWeights)))
    ifelse(dim(summ)[1] > 1, summ[2, 4], 1)
  }, Y = Y, family = family)
  ## rank the p-values; give age, BMI, bhvrisk the lowest rank (will always set to TRUE anyways)
  ranked_vars <- rank(listp, ties = "average")
  ranked_vars[names(X) %in% c("age", "BMI", "bhvrisk")] <- 999
  return(ranked_vars)
}

# Screening 1. Dynamic range #
screen_dynamic_range_plus_exposure <- function(Y, X, family, obsWeights, nVar = 4, ...) {
  # set all vars to false
  vars <- rep(FALSE, ncol(X))
  
  # keep only those with dynamic range: 20th percentile != 80th percentile
  x_quantiles <- apply(X, 2, function(x) quantile(x, probs = c(0.2, 0.8)))
  vars <- apply(x_quantiles, 2, function(x) round(x[1], 4) != round(x[2], 4))
  # also keep the first three columns of X (correspond to age, BMI, bhvrisk)
  vars[names(X) %in% c("age", "BMI", "bhvrisk")] <- TRUE
  # keep only a max of nVar immune markers; rank by univariate p-value in a model adjusting for age, BMI, bhvrisk
  X_initial_screen <- X %>% select(names(X)[vars], "age", "BMI", "bhvrisk")
  ranked_vars <- rank_univariate_logistic_pval_plus_exposure(Y, X_initial_screen, family, obsWeights )
  vars[vars][ranked_vars > nVar] <- FALSE
  
  # also keep the first three columns of X (correspond to age, BMI, bhvrisk)
  vars[names(X) %in% c("age", "BMI", "bhvrisk")] <- TRUE
  return(vars)
}

# Screening 2. Dynamic range score #
screen_dynamic_range_score_plus_exposure <- function(Y, X, family, obsWeights, var_super = var.super, nVar = 4, ...) {
  # set all to false
  vars <- rep(FALSE, ncol(X))
  # need to apply with the correct label in place of X
  vars_sd_ratio <- ifelse(is.na(var_super$sd.ratio), TRUE, var_super$sd.ratio > quantile(var_super$sd.ratio, probs = c(0.5), na.rm = TRUE))
  vars <- names(X) %in% var_super$varname[vars_sd_ratio] | names(X) %in% paste0(var_super$varname[vars_sd_ratio], "_bin")
  names(vars) <- colnames(X)
  # also keep the first three columns of X (correspond to age, BMI, bhvrisk)
  vars[names(X) %in% c("age", "BMI", "bhvrisk")] <- TRUE
  # keep only a max of nVar immune markers; rank by univariate p-value in a model adjusting for age, BMI, bhvrisk
  X_initial_screen <- X %>% select(names(X)[vars], "age", "BMI", "bhvrisk")
  ranked_vars <- rank_univariate_logistic_pval_plus_exposure(Y, X_initial_screen, family, obsWeights)
  vars[vars][ranked_vars > nVar] <- FALSE
  # also keep the first three columns of X (correspond to age, BMI, bhvrisk)
  vars[names(X) %in% c("age", "BMI", "bhvrisk")] <- TRUE
  return(vars)
}

# Screening 3. Lasso #
screen_lasso_plus_exposure <- function(Y, X, family, obsWeights, alpha = 1, nVar = 500, ...) {
  # Lasso #
  set.seed(123)
  #res.ls <- cv.glmnet( x = as.matrix(X), y =  as.matrix(Y), family = "binomial" , type.measure = 'auc', nfolds = 5, alpha = 1 ) # Lasso penalty
  res.ls <- cv.glmnet( x = as.matrix(X), y =  as.matrix(Y), weights = obsWeights, family = "binomial" , type.measure = 'auc', nfolds = 5, alpha = 1 ) # Lasso penalty
  vars.ls <- (coef( res.ls, s = res.ls$lambda.min ) != 0)[-1]
  vars <- vars.ls
  names(vars) <- colnames(X)
  
  # also keep the first three columns of X (correspond to age, BMI, bhvrisk)
  vars[names(X) %in% c("age", "BMI", "bhvrisk")] <- TRUE
  # keep only a max of nVar immune markers; rank by univariate p-value in a model adjusting for age, BMI, bhvrisk
  X_initial_screen <- X %>% select(names(X)[vars], "age", "BMI", "bhvrisk")
  ranked_vars <- rank_univariate_logistic_pval_plus_exposure(Y, X_initial_screen, family, obsWeights)
  vars[vars][ranked_vars > nVar] <- FALSE
  # also keep the first three columns of X (correspond to age, BMI, bhvrisk)
  vars[names(X) %in% c("age", "BMI", "bhvrisk")] <- TRUE
  return(vars)
}

# Screening 4. High correlation #
screen_highcor_plus_exposure <- function(Y, X, family, obsWeights, nVar = 4, ...) {
  # set all vars to FALSE
  vars <- rep(FALSE, ncol(X))
  # compute pairwise correlations between all marker vars
  cors <- cor(X, method = "spearman")
  diag(cors) <- NA
  cor_less_0.9 <- (cors <= 0.9)
  # screen out those with r > 0.9
  vars <- apply(cor_less_0.9, 1, function(x) all(x, na.rm = TRUE))
  # also keep the first three columns of X (correspond to age, BMI, bhvrisk)
  vars[names(X) %in% c("age", "BMI", "bhvrisk")] <- TRUE
  # keep only a max of nVar immune markers; rank by univariate p-value in a model adjusting for age, BMI, bhvrisk
  X_initial_screen <- X %>% select(names(X)[vars], "age", "BMI", "bhvrisk")
  ranked_vars <- rank_univariate_logistic_pval_plus_exposure(Y, X_initial_screen, family, obsWeights)
  vars[vars][ranked_vars > nVar] <- FALSE
  ## make sure that age, BMI, bhvrisk are true
  vars[names(X) %in% c("age", "BMI", "bhvrisk")] <- TRUE
  return(vars)
}

# Screening 5. Univariate logistic regression (ULR) #
screen_ulr_pval_plus_exposure <- function(Y, X, family, obsWeights, minPvalue=0.05, minscreen = 2, nVar = 4, ...) {
  ## logistic regression of outcome on each variable
  listp <- apply(X, 2, function(x, Y, family) {
    summ <- coef(summary(glm(Y ~ x + X$age + X$BMI + X$bhvrisk, family = family, weights = obsWeights)))
    ifelse(dim(summ)[1] > 1, summ[2, 4], 1)
  }, Y = Y, family = family)
  vars.ulr <- (listp <= minPvalue)
  
  # also keep the first three columns of X (correspond to age, BMI, bhvrisk)
  vars <- vars.ulr
  vars[names(X) %in% c("age", "BMI", "bhvrisk")] <- TRUE
  if (sum(vars) < minscreen) {
    warning("number of variables with p value less than minPvalue is less than minscreen")
    vars[rank(listp) <= minscreen] <- TRUE
  }
  # keep only a max of nVar immune markers; rank by univariate p-value in a model adjusting for age, BMI, bhvrisk
  X_initial_screen <- X %>% select(names(X)[vars], "age", "BMI", "bhvrisk")
  ranked_vars <- rank(listp[vars], ties = "average")
  ranked_vars[names(ranked_vars) %in% c("age", "BMI", "bhvrisk")] <- 999
  # ranked_vars <- rank_univariate_logistic_pval_plus_exposure(Y, X_initial_screen, family, obsWeights)
  vars[vars][ranked_vars > nVar] <- FALSE
  vars[names(X) %in% c("age", "BMI", "bhvrisk")] <- TRUE
  return(vars)
}

# Total screening function #
screen.dat.index <- function(Y, X, X_markers, fit.set = c('None','BAMA','Tcells','BAMA_Tcells','All'),
                             screen.dat.method = c('all','dr','drs','ls','hc','ulr'), 
                             screen.index.method = c('all','dr','drs','ls','hc','ulr'), obsWeights){
  
  # Candidate set #
  if( fit.set == 'None' ){
    # 2-1. No markers (only clinical covariates : age, BMI, bhrisk) #
    var_set_none <- rep(FALSE, ncol(X_markers))
    var_set_none <- c( rep(TRUE, 3), var_set_none )
    dat <- cbind(Y = Y, X[,var_set_none])
  } else if( fit.set == 'BAMA' ){
    # 2-2. BAMA markers (IgG + IgA + IgG3) #
    var_set_igg_iga_igg3 <- get_nms_group_all_antigens(X_markers, assays = c("IgG", "IgA", "IgG3"))
    var_set_igg_iga_igg3 <- c( rep(TRUE, 3), var_set_igg_iga_igg3 )
    dat <- cbind(Y = Y, X[,var_set_igg_iga_igg3])
  } else if( fit.set == 'Tcells' ){
    # 2-3. T cells (CD4 and CD8) #
    var_set_tcells <- get_nms_group_all_antigens(X_markers, assays = c("CD4", "CD8"))
    var_set_tcells <- c( rep(TRUE, 3), var_set_tcells )
    dat <- cbind(Y = Y, X[,var_set_tcells])
  } else if( fit.set == 'BAMA_Tcells' ){
    # 2-5. BAMA markers + T cells #
    var_set_igg_iga_igg3_tcells <- get_nms_group_all_antigens(X_markers, assays = c("IgG", "IgA", "IgG3", "CD4", "CD8"))
    var_set_igg_iga_igg3_tcells <- c( rep(TRUE, 3), var_set_igg_iga_igg3_tcells )
    dat <- cbind(Y = Y, X[,var_set_igg_iga_igg3_tcells])
  } else if( fit.set == 'All' ){
    # 2-4. All markers #
    var_set_all <- rep(TRUE, ncol(X_markers))
    var_set_all <- c( rep(TRUE, 3), var_set_all )
    dat <- cbind(Y = Y, X[,var_set_all])
  }
  
  # Screened data #
  if( screen.dat.method == 'all' ){
    screen.dat.var <- rep(TRUE, (ncol(dat)-1))
  } else if( screen.dat.method == 'dr' ){
    screen.dat.var <- screen_dynamic_range_plus_exposure( Y = dat$Y, X = dat[,colnames(dat)!='Y'], family = 'binomial', obsWeights = obsWeights )
  } else if( screen.dat.method == 'drs' ){
    screen.dat.var <- screen_dynamic_range_score_plus_exposure( Y = dat$Y, X = dat[,colnames(dat)!='Y'], family = 'binomial', obsWeights = obsWeights )
  } else if( screen.dat.method == 'ls' ){
    screen.dat.var <- screen_lasso_plus_exposure( Y = dat$Y, X = dat[,colnames(dat)!='Y'], family = 'binomial', obsWeights = obsWeights )
  } else if( screen.dat.method == 'hc' ){
    screen.dat.var <- screen_highcor_plus_exposure( Y = dat$Y, X = dat[,colnames(dat)!='Y'], family = 'binomial', obsWeights = obsWeights )
  } else if( screen.dat.method == 'ulr' ){
    screen.dat.var <- screen_ulr_pval_plus_exposure( Y = dat$Y, X = dat[,colnames(dat)!='Y'], family = 'binomial', obsWeights = obsWeights )
  }
  
  # Screened index #
  if( screen.index.method == 'all' ){
    screen.index.var <- rep(TRUE, (ncol(dat)-1))
  } else if( screen.index.method == 'dr' ){
    screen.index.var <- screen_dynamic_range_plus_exposure( Y = dat$Y, X = dat[,colnames(dat)!='Y'], family = 'binomial', obsWeights = obsWeights )
  } else if( screen.index.method == 'drs' ){
    screen.index.var <- screen_dynamic_range_score_plus_exposure( Y = dat$Y, X = dat[,colnames(dat)!='Y'], family = 'binomial', obsWeights = obsWeights )
  } else if( screen.index.method == 'ls' ){
    screen.index.var <- screen_lasso_plus_exposure( Y = dat$Y, X = dat[,colnames(dat)!='Y'], family = 'binomial', obsWeights = obsWeights )
  } else if( screen.index.method == 'hc' ){
    screen.index.var <- screen_highcor_plus_exposure( Y = dat$Y, X = dat[,colnames(dat)!='Y'], family = 'binomial', obsWeights = obsWeights )
  } else if( screen.index.method == 'ulr' ){
    screen.index.var <- screen_ulr_pval_plus_exposure( Y = dat$Y, X = dat[,colnames(dat)!='Y'], family = 'binomial', obsWeights = obsWeights )
  }
  
  dat.X <- list( case = dat[dat$Y == 1, c(FALSE, screen.dat.var)], control = dat[dat$Y == 0, c(FALSE, screen.dat.var)])
  return(list(screen.index.var = screen.index.var, dat = dat.X))
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

# RF-based CV-AUC #
get.rf.cvauc = function(dat, cv.scheme, obsWeights, method=c('sRF','sRF_under','sRF_over','tRF'), seed=1){
  
  # Split data for k-fold CV #
  splits <- get.splits(dat, cv.scheme, seed)
  
  # Standard Random Forest (sRF) #
  if( method == 'sRF' ){
    cv.aucs <-  mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      weights.train <- obsWeights[as.numeric(rownames(dat.train))] ; weights.test <- obsWeights[as.numeric(rownames(dat.test))]
      
      set.seed( 123 )
      fit.rf <- ranger( factor(Y)~., data = dat.train, case.weights = weights.train, probability = TRUE, min.node.size = 1 )
      pred.rf <- predict( fit.rf, data = dat.test )
      measure_auc( pred.rf$predictions[,'1'], dat.test$Y, weights = weights.test )$point_est
    }, mc.cores = 4 )
    cv.aucs <- unlist( cv.aucs )
  }
  
  # sRF + under-sampling (sRF_under) #
  if( method == 'sRF_under' ){
    cv.aucs <-  mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      weights.train <- obsWeights[as.numeric(rownames(dat.train))] ; weights.test <- obsWeights[as.numeric(rownames(dat.test))]
      weights.sampling <- rep(NA, nrow(dat.train))
      weights.sampling[dat.train$Y == 1] <- 5 ; weights.sampling[dat.train$Y == 0] <- 1
      
      set.seed( 123 )
      fit.rf <- ranger( factor(Y)~., data = dat.train, probability = TRUE, replace = TRUE, case.weights = weights.sampling, 
                        sample.fraction = 40/120, keep.inbag = TRUE, min.node.size = 1 )
      pred.rf <- predict( fit.rf, data = dat.test)
      measure_auc( pred.rf$predictions[,'1'], dat.test$Y, weights = weights.test )$point_est
    }, mc.cores = 4 )
    cv.aucs <- unlist( cv.aucs )
  }
  
  # sRF + over-sampling (sRF_over) #
  if( method == 'sRF_over' ){
    cv.aucs <-  mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      set.seed( 123 )
      idx1 <- sample( 1:sum(dat.train$Y==1), sum(dat.train$Y==0), replace = TRUE ) # Oversampling minority class
      dat.temp <- list( case = subset(dat.train, Y==1, select=-Y), control = subset(dat.train, Y==0, select=-Y) )
      dat.train.mod <- rbind( data.frame(Y=1, dat.temp$case[idx1,,drop=F]),   data.frame(Y=0, dat.temp$control) )
      weights.train <- obsWeights[floor(as.numeric(rownames(dat.train.mod)))] ; weights.test <- obsWeights[as.numeric(rownames(dat.test))]
      
      set.seed( 123 )
      fit.rf <- ranger( factor(Y)~., data = dat.train.mod, probability = TRUE, replace = TRUE, case.weights = weights.train, min.node.size = 1 )
      pred.rf <- predict( fit.rf, data = dat.test)
      measure_auc( pred.rf$predictions[,'1'], dat.test$Y, weights = weights.test )$point_est
    }, mc.cores = 4 )
    cv.aucs <- unlist( cv.aucs )
  }
  
  # Tuned Random Forest (tRF) #
  if( method == 'tRF' ){
    cv.aucs <-  mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      dat.train$Y <- as.factor(dat.train$Y) ; dat.test$Y <- as.factor(dat.test$Y)
      weights.train <- obsWeights[as.numeric(rownames(dat.train))] ; weights.test <- obsWeights[as.numeric(rownames(dat.test))]
      
      set.seed( 123 )
      rf.task <- makeClassifTask(data = dat.train, target = 'Y')
      res.tunerf <- tuneRanger( rf.task, measure = list(auc), num.trees = 500, iters.warmup = 50, iters = 100, save.file.path = NULL, 
                                tune.parameters = c("mtry", "min.node.size","sample.fraction"), parameters = list(replace = TRUE) )
      pred.tunerf <- predict( res.tunerf$model, newdata = dat.test )
      measure_auc( pred.tunerf$data$prob.1, dat.test$Y, weights = weights.test )$point_est
    }, mc.cores = 4 )
    cv.aucs <- unlist( cv.aucs )
  }
  
  cv.aucs
}

# Stacking CV-AUC #
get.st.cvauc = function(dat, cv.scheme, obsWeights, var.index, method, seed=1){
  
  # Outer layer: 5-fold CV #
  splits <- get.splits(dat, cv.scheme, seed)
  
  # Inner layer: 5-fold CV #
  my_control <- trainControl(
    method="cv",
    number=119,
    savePredictions="final",
    classProbs=TRUE,
    summaryFunction=twoClassSummary
  )
  
  # Method of estimating regression coefficients for stacking #
  method <- get(method, mode = 'function')()
  if(!is.null(method$require)) {
    sapply(method$require, function(x) require(force(x), character.only = TRUE))
  }
  
  # Fit candidate learners #
  cv.aucs <-  mclapply( splits, function(split){
    # The target variables must be named as 'Y' #
    dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
    dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
    # Y=1, 0 must be converted to 'case' and 'contrl' factor levels respectively #
    dat.train$Y <- as.factor(dat.train$Y) ; levels(dat.train$Y) <- c('control','case')
    weights.train <- obsWeights[as.numeric(rownames(dat.train))] ; weights.test <- obsWeights[as.numeric(rownames(dat.test))]
    
    # Stacking with two candidate learners #
    if(length(var.index) == 2){
      # For sRF #
      rf_grid <- expand.grid(mtry = floor(sqrt(sum(var.index$rf)-1)), splitrule = 'gini', min.node.size = 1)
      set.seed( 123 )
      model_list <- caretList(
        list(Y~., data=dat.train[,var.index$glm], weights = weights.train), 
        list(Y~., data=dat.train[,var.index$rf], tuneGrid = rf_grid, weights = weights.train),
        trControl=my_control,
        methodList=c('glm', 'ranger')
      )
    }
    
    # Stacking with three candidate learners #
    if(length(var.index) == 3){
      # For sRF #
      rf_grid <- expand.grid(mtry = floor(sqrt(sum(var.index$rf)-1)), splitrule = 'gini', min.node.size = 1)
      rf.1_grid <- expand.grid(mtry = floor(sqrt(sum(var.index$rf.1)-1)), splitrule = 'gini', min.node.size = 1)
      set.seed( 123 )
      model_list <- caretList(
        list(Y~., data=dat.train[,var.index$glm], weights = weights.train), 
        list(Y~., data=dat.train[,var.index$rf], tuneGrid = rf_grid, weights = weights.train),
        list(Y~., data=dat.train[,var.index$rf.1], tuneGrid = rf.1_grid, weights = weights.train),
        trControl=my_control,
        methodList=c('glm', 'ranger', 'ranger')
      )
    }
    
    # Super learning #
    set.seed( 123 )
    pred.cv <- makePredObsMatrix(model_list)
    res.st.fit <- method$computeCoef(Z = pred.cv$preds, Y =(as.numeric(pred.cv$obs)-1), obsWeights = weights.train, 
                                     libraryNames = names(var.index), verbose = FALSE)
    pred.test <- predict.caretList(model_list, newdata=dat.test)
    pred.st <- method$computePred(predY = pred.test, coef = res.st.fit$coef)
  
    # Predict stacking on test set #
    measure_auc( pred.st, dat.test$Y, weights = weights.test )$point_est
  }, mc.cores = 4 )
  cv.aucs <- unlist( cv.aucs )
  
  cv.aucs
}


#############################
# more general functions

screen_lasso <- function(Y, X, family, obsWeights=rep(1, nrow(X)), alpha = 1) {
  set.seed(123)
  res.ls <- cv.glmnet( x = as.matrix(X), y =  as.matrix(Y), weights = obsWeights, family = "binomial" , type.measure = 'auc', nfolds = 5, alpha = alpha ) # Lasso penalty
  vars.ls <- (coef( res.ls, s = res.ls$lambda.min ) != 0)[-1]
  vars <- vars.ls
  names(vars) <- colnames(X)
  return(vars)
}

screen_ulr <- function(Y, X, family, obsWeights=rep(1, nrow(X)), cutoff=0.1) {
  pvals=sapply (1:ncol(X), function(i) {
    fit=glm.fit(cbind(1, X[,i]), Y, family=binomial(), weights=obsWeights)
    pval=last(c(summary.glm(fit)$coef)) # glm.fit is faster than glm
    pval
  })
  sel=pvals<cutoff
  sel
}
