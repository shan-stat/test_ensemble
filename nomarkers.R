# 1. Re-define function for only using baseline variables #
no.markers.cvauc = function(dat, cv.scheme, obsWeights, measure=c('trCV-AUC','tCV-AUC'), method=c('nodesize','mtry','sampsize'), value, seed=1){
  
  # Split data for k-fold CV #
  splits <- get.splits(dat, cv.scheme, seed)
  
  # nodesize #
  if( measure == 'trCV-AUC' & method == 'nodesize' ){
    cv.aucs <-  mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      weights.train <- obsWeights[as.numeric(rownames(dat.train))] ; weights.test <- obsWeights[as.numeric(rownames(dat.test))]
    
      # Training AUC #
      set.seed( 123 )
      fit.rf <- ranger( factor(Y)~., data = dat.train, case.weights = weights.train, probability = TRUE, min.node.size = value )
      pred.rf <- predict( fit.rf, data = dat.train )
      measure_auc( pred.rf$predictions[,'1'], dat.train$Y, weights = weights.train )$point_est
    }, mc.cores = 4 )
    cv.aucs <- unlist( cv.aucs )
  }
  
  if( measure == 'tCV-AUC' & method == 'nodesize' ){
    cv.aucs <-  mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      weights.train <- obsWeights[as.numeric(rownames(dat.train))] ; weights.test <- obsWeights[as.numeric(rownames(dat.test))]
      
      # Test AUC #
      set.seed( 123 )
      fit.rf <- ranger( factor(Y)~., data = dat.train, case.weights = weights.train, probability = TRUE,  min.node.size = value )
      pred.rf <- predict( fit.rf, data = dat.test )
      measure_auc( pred.rf$predictions[,'1'], dat.test$Y, weights = weights.test )$point_est
    }, mc.cores = 4 )
    cv.aucs <- unlist( cv.aucs )
  }
  
  # mtry #
  if( measure == 'trCV-AUC' & method == 'mtry' ){
    cv.aucs <-  mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      weights.train <- obsWeights[as.numeric(rownames(dat.train))] ; weights.test <- obsWeights[as.numeric(rownames(dat.test))]
      
      # Training AUC #
      set.seed( 123 )
      fit.rf <- ranger( factor(Y)~., data = dat.train, case.weights = weights.train, probability = TRUE, min.node.size = 1, mtry = value )
      pred.rf <- predict( fit.rf, data = dat.train)
      measure_auc( pred.rf$predictions[,'1'], dat.train$Y, weights = weights.train )$point_est
    }, mc.cores = 4 )
    cv.aucs <- unlist( cv.aucs )
  }

  if( measure == 'tCV-AUC' & method == 'mtry' ){
    cv.aucs <-  mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      weights.train <- obsWeights[as.numeric(rownames(dat.train))] ; weights.test <- obsWeights[as.numeric(rownames(dat.test))]
      
      # Test AUC #
      set.seed( 123 )
      fit.rf <- ranger( factor(Y)~., data = dat.train, case.weights = weights.train, probability = TRUE, min.node.size = 1, mtry = value )
      pred.rf <- predict( fit.rf, data = dat.test )
      measure_auc( pred.rf$predictions[,'1'], dat.test$Y, weights = weights.test )$point_est
    }, mc.cores = 4 )
    cv.aucs <- unlist( cv.aucs )
  }
  
  # sampsize #
  if( measure == 'trCV-AUC' & method == 'sampsize' ){
    cv.aucs <-  mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      weights.train <- obsWeights[as.numeric(rownames(dat.train))] ; weights.test <- obsWeights[as.numeric(rownames(dat.test))]
      
      # Training AUC #
      set.seed( 123 )
      fit.rf <- ranger( factor(Y)~., data = dat.train, case.weights = weights.train, probability = TRUE, 
                        min.node.size = 1, sample.fraction = value )
      pred.rf <- predict( fit.rf, data = dat.train )
      measure_auc( pred.rf$predictions[,'1'], dat.train$Y, weights = weights.train )$point_est
    }, mc.cores = 4 )
    cv.aucs <- unlist( cv.aucs )
  }
  
  if( measure == 'tCV-AUC' & method == 'sampsize' ){
    cv.aucs <-  mclapply( splits, function(split){
      dat.train <- rbind( data.frame(Y=1, dat$case[split$training$case,,drop=F]),   data.frame(Y=0, dat$control[split$training$control,,drop=F]) )
      dat.test <- rbind( data.frame(Y=1, dat$case[split$test$case,,drop=F]),       data.frame(Y=0, dat$control[split$test$control,,drop=F]) )
      weights.train <- obsWeights[as.numeric(rownames(dat.train))] ; weights.test <- obsWeights[as.numeric(rownames(dat.test))]
      
      # Test AUC #
      set.seed( 123 )
      fit.rf <- ranger( factor(Y)~., data = dat.train, case.weights = weights.train, probability = TRUE, 
                        min.node.size = 1, sample.fraction = value )
      pred.rf <- predict( fit.rf, data = dat.test )
      measure_auc( pred.rf$predictions[,'1'], dat.test$Y, weights = weights.test )$point_est
    }, mc.cores = 4 )
    cv.aucs <- unlist( cv.aucs )
  }
  
  cv.aucs
}



# 1. Data setting #
data( 'dat.505', package = 'HVTN505' )
suppressWarnings( data( 'var.super', package = 'HVTN505' ) )

# Pre-scaling each variables: mean 0 and sd 1 #
for( i in var.super$varname ){
  dat.505[[i]] <- scale( dat.505[[i]], center = mean(dat.505[[i]][dat.505$trt == 1]), scale = sd(dat.505[[i]][dat.505$trt == 1]))
  dat.505[[i %.% '_bin']] <- scale( dat.505[[i %.% '_bin']], center = mean(dat.505[[i %.% '_bin']][dat.505$trt == 1]), scale = sd(dat.505[[i %.% '_bin']][dat.505$trt == 1]))
}
for( i in c( 'age', 'BMI', 'bhvrisk' ) ) {
  dat.505[[i]] <- scale( dat.505[[i]], center = mean(dat.505[[i]][dat.505$trt == 1]), scale = sd(dat.505[[i]][dat.505$trt == 1]))
}

# Raw data (vaccinees + placebo recipients = 189 obs) #
X_markers <- dat.505 %>% select( var.super$varname, paste0(var.super$varname, '_bin') )
X_covariates <- dat.505 %>% select( age, BMI, bhvrisk )
X <- data.frame( trt = dat.505$trt, X_covariates, X_markers )
X_full_none <- data.frame( trt = dat.505$trt, X_covariates )
weights <- dat.505$wt
Y <- dat.505$case

# Vaccinees (150 obs) #
vaccinees <- cbind.data.frame( Y, weights, X ) %>% filter( trt == 1 ) %>% select( -trt )
X_none <- cbind.data.frame( Y, weights, X_full_none ) %>% filter( trt == 1 ) %>% select( -trt, -Y, -weights )
Y_vaccine <- vaccinees$Y
weights_vaccine <- vaccinees$weights
X_vaccine <- vaccinees %>% select( -Y, -weights )

# None (only covariates : age, BMI, bhrisk) #
assays <- unique(var.super$assay)
antigens <- unique(var.super$antigen)

var_set_none <- rep(FALSE, ncol(X_markers))
var_set_none <- c( rep(TRUE, 3), var_set_none )
dat <- cbind( Y = Y_vaccine, X_vaccine[,var_set_none] )
dat.X <- list( case = dat[dat$Y == 1, -1], control = dat[dat$Y == 0, -1])



# 2. Fitting codes #
# 2-1. Mtry #
mtry.mat <- matrix( NA, nrow = 10, ncol = 3 )
colnames( mtry.mat ) <- paste0('mtry_', 1:3)
mtry.can <- 1:3

for( i in 1:length(mtry.can) ){
  print(i)
  for( j in 1:10 ){
    print(j)
    seed <- j
    res.temp <- no.markers.cvauc( dat = dat.X, cv.scheme = '5fold', obsWeights = weights_vaccine, 
                                  measure = 'trCV-AUC', method = 'mtry', value = mtry.can[i], seed = seed ) 
    mtry.mat[j, i] <- mean( res.temp )
  }
}
apply(mtry.mat, 2, mean)
write.table( mtry.mat, '/Users/shan/Desktop/Paper/YFong/2.HVTN/Result/mtry_tauc.txt', col.names = TRUE, sep =',' )



# 2-2. Nodesize #
ndsize.mat <- matrix( NA, nrow = 10, ncol = 11 )
colnames( ndsize.mat ) <- c( 'ndsize_1', paste0('ndsize_', seq(0.1, 0.9, 0.1), 'n'), 'ndsize_1.0n')
nodesize.can <- c(1, seq(0.1, 1, 0.1)*120 )

for( i in 1:length(nodesize.can) ){
  print(i)
  for( j in 1:10 ){
    print(j)
    seed <- j
    res.temp <- no.markers.cvauc( dat = dat.X, cv.scheme = '5fold', obsWeights = weights_vaccine, 
                                  measure = 'tCV-AUC', method = 'nodesize', value = nodesize.can[i], seed = seed ) 
    ndsize.mat[j, i] <- mean( res.temp )
  }
}
apply(ndsize.mat, 2, mean)
write.table( ndsize.mat, '/Users/shan/Desktop/Paper/YFong/2.HVTN/Result/nodesize_tauc.txt', col.names = TRUE, sep =',' )



# 2-3. Sampsize #
spsize.mat <- matrix( NA, nrow = 10, ncol = 10 )
colnames( spsize.mat ) <- c( paste0('spsize_', seq(0.1, 0.9, 0.1), 'n'), 'spsize_1.0n')
sampsize.can <- (seq(0.1, 1, 0.1)*120)/120

for( i in 1:length(sampsize.can) ){
  print(i)
  for( j in 1:10 ){
    print(j)
    seed <- j
    res.temp <- no.markers.cvauc( dat = dat.X, cv.scheme = '5fold', obsWeights = weights_vaccine, 
                                  measure = 'tCV-AUC', method = 'sampsize', value = sampsize.can[i], seed = seed ) 
    spsize.mat[j, i] <- mean( res.temp )
  }
}
apply(spsize.mat, 2, mean)
write.table( spsize.mat, '/Users/shan/Desktop/Paper/YFong/2.HVTN/Result/spsize_tauc.txt', col.names = TRUE, sep =',' )



# 3. Summary results #
train.auc <- read.table( '/Users/shan/Desktop/Paper/YFong/2.HVTN/Result/sRF_cvauc.txt', header = T, sep = ',' )
round(apply( train.auc, 2, mean ), 3)
