# Generating data #
data( 'dat.505', package = 'HVTN505' )
suppressWarnings( data( 'var.super', package = 'HVTN505' ) )

# Scaling vaccine recipients to have mean 0, sd 1 for all vars #
for( i in var.super$varname ){
  dat.505[[i]] <- scale( dat.505[[i]], center = mean(dat.505[[i]][dat.505$trt == 1]), scale = sd(dat.505[[i]][dat.505$trt == 1]))
  dat.505[[i %.% '_bin']] <- scale( dat.505[[i %.% '_bin']], center = mean(dat.505[[i %.% '_bin']][dat.505$trt == 1]), scale = sd(dat.505[[i %.% '_bin']][dat.505$trt == 1]))
}
for( i in c( 'age', 'BMI', 'bhvrisk' ) ) {
  dat.505[[i]] <- scale( dat.505[[i]], center = mean(dat.505[[i]][dat.505$trt == 1]), scale = sd(dat.505[[i]][dat.505$trt == 1]))
}

X_markers <- dat.505 %>% select( var.super$varname, paste0(var.super$varname, '_bin') )
X_covariates <- dat.505 %>% select( age, BMI, bhvrisk )
X <- data.frame( trt = dat.505$trt, X_covariates, X_markers )
X_full_none <- data.frame( trt = dat.505$trt, X_covariates )
weights <- dat.505$wt
Y <- dat.505$case

vaccinees <- cbind.data.frame( Y, weights, X ) %>% filter( trt == 1 ) %>% select( -trt )
X_none <- cbind.data.frame( Y, weights, X_full_none ) %>% filter( trt == 1 ) %>% select( -trt, -Y, -weights )
Y_vaccine <- vaccinees$Y
weights_vaccine <- vaccinees$weights
X_vaccine <- vaccinees %>% select( -Y, -weights )



### Fitting RF ###
assays <- unique(var.super$assay)
antigens <- unique(var.super$antigen)

pred.mat <- matrix( NA, nrow = 10, ncol = 11 )
colnames( pred.mat ) <- c( 'None', 'IgG_IgA', 'IgG3', 'Tcells', 'FxAb', 'IgG_IgA_IgG3', 'IgG_IgA_Tcells',
                           'IgG_IgA_IgG3_Tcells', 'IgG_IgA_IgG3_FxAb', 'Tcells_FxAb', 'All' )


## 1. None (only covariates : age, BMI, bhrisk) ##
var_set_none <- rep(FALSE, ncol(X_markers))
var_set_none <- c( rep(TRUE, 3), var_set_none )
dat <- cbind( Y = Y_vaccine, X_vaccine[,var_set_none] )
dat.X <- list( case = dat[dat$Y == 1, -1], control = dat[dat$Y == 0, -1])

# K-fold CV #
for( i in 1:10 ){
  print(i)
  seed <- i
  res.temp <- get.kfold.cvauc( dat = dat.X, cv.scheme = '5fold', method = 'RF', seed = seed ) # 5-fold
  pred.mat[i, 'None'] <- mean( res.temp )
}
pred.mat
