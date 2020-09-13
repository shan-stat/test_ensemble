### Required file for Stacking ###
setwd('/Users/shan/Desktop/caretEnsemble_mod/R')
source('helper_functions.R')
source('caretList.R')
source('caretEnsemble.R')
source('caretStack.R')



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



# 2. Screening #
assays <- unique(var.super$assay)
antigens <- unique(var.super$antigen)

pred.mat <- matrix( NA, nrow = 10, ncol = 5 )
colnames( pred.mat ) <- c( 'None', 'BAMA', 'Tcells', 'All', 'BAMA_Tcells' )

# 2-1. No markers (only clinical covariates : age, BMI, bhrisk) #
var_set_none <- rep(FALSE, ncol(X_markers))
var_set_none <- c( rep(TRUE, 3), var_set_none )
dat <- cbind( Y = Y_vaccine, X_vaccine[,var_set_none] )
dat.X <- screen.index(dat = dat, method = 'all', obsWeights = weights_vaccine)$dat
screen.var <- screen.index(dat = dat, method = 'all', obsWeights = weights_vaccine)$screen.var

# 2-2. BAMA markers (IgG + IgA + IgG3) #
var_set_igg_iga_igg3 <- get_nms_group_all_antigens(X_markers, assays = c("IgG", "IgA", "IgG3"))
var_set_igg_iga_igg3 <- c( rep(TRUE, 3), var_set_igg_iga_igg3 )
dat <- cbind( Y = Y_vaccine, X_vaccine[,var_set_igg_iga_igg3] )
dat.X <- screen.index(dat = dat, method = 'all', obsWeights = weights_vaccine)$dat
screen.var <- screen.index(dat = dat, method = 'ls', obsWeights = weights_vaccine)$screen.var

# 2-3. T cells (CD4 and CD8) #
var_set_tcells <- get_nms_group_all_antigens(X_markers, assays = c("CD4", "CD8"))
var_set_tcells <- c( rep(TRUE, 3), var_set_tcells )
dat <- cbind( Y = Y_vaccine, X_vaccine[,var_set_tcells] )
dat.X <- screen.index(dat = dat, method = 'all', obsWeights = weights_vaccine)$dat
screen.var <- screen.index(dat = dat, method = 'ls', obsWeights = weights_vaccine)$screen.var

# 2-4. All markers #
var_set_all <- rep(TRUE, ncol(X_markers))
var_set_all <- c( rep(TRUE, 3), var_set_all )
dat <- cbind( Y = Y_vaccine, X_vaccine[,var_set_all] )
dat.X <- screen.index(dat = dat, method = 'all', obsWeights = weights_vaccine)$dat
screen.var <- screen.index(dat = dat, method = 'ls', obsWeights = weights_vaccine)$screen.var

# 2-5. BAMA markers + T cells #
var_set_igg_iga_igg3_tcells <- get_nms_group_all_antigens(X_markers, assays = c("IgG", "IgA", "IgG3", "CD4", "CD8"))
var_set_igg_iga_igg3_tcells <- c( rep(TRUE, 3), var_set_igg_iga_igg3_tcells )
dat <- cbind( Y = Y_vaccine, X_vaccine[,var_set_igg_iga_igg3_tcells] )
dat.X <- screen.index(dat = dat, method = 'all', obsWeights = weights_vaccine)$dat
screen.var <- screen.index(dat = dat, method = 'ls', obsWeights = weights_vaccine)$screen.var



# 3. Fitting codes #
for( i in 1:10 ){
  print(i)
  seed <- i
  
  # RF-based models #
  res.temp <- get.rf.cvauc( dat = dat.X, cv.scheme = '5fold', obsWeights = weights_vaccine, method = 'sRF', seed = seed )
  
  # Stacking with 3 learners #
  #res.temp <- get.st.cvauc( dat = dat.X, cv.scheme = '5fold', obsWeights = weights_vaccine, 
  #                          var.index = list(glm  = c(TRUE, rep(TRUE,3), rep(FALSE,ncol(dat.X$case)-3)), # clinical covariates
  #                                           rf   = c(TRUE, screen.var),                                  # clinical covariates + screened variables
  #                                           rf.1 = c(TRUE, rep(TRUE,ncol(dat.X$case)))),               # all variables
  #                          seed = seed )
  
  # Stacking with 2 learners #
  #res.temp <- get.st.cvauc( dat = dat.X, cv.scheme = '5fold', obsWeights = weights_vaccine, 
  #                          var.index = list(glm = c(TRUE, rep(TRUE,3), rep(FALSE,ncol(dat.X$case)-3)), # clinical covariates
  #                                           rf  = c(TRUE, screen.var)),                                  # clinical covariates + screened variables
  #                          seed = seed )
  
  pred.mat[i, 'None']        <- mean( res.temp )
  #pred.mat[i, 'BAMA']        <- mean( res.temp )
  #pred.mat[i, 'Tcells']      <- mean( res.temp )
  #pred.mat[i, 'All']         <- mean( res.temp )
  #pred.mat[i, 'BAMA_Tcells'] <- mean( res.temp )
}
apply(pred.mat, 2, mean)

write.table( pred.mat, file = '/Users/shan/Desktop/Paper/YFong/2.HVTN/Result/ST_can_2.txt', col.names = T, sep = ',' )
