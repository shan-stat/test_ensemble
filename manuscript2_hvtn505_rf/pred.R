source('helper_cvauc.R')
source('helper_rf_hvtn.R')

### For Stacking ###
setwd('/Users/shan/Desktop/Paper/YFong/2.HVTN/Code/caretEnsemble_modified/R')
source('helper_functions.R')
source('caretList.R')
source('caretEnsemble.R')
source('caretStack.R')
source('method.R')



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

# for RF, screen.dat.method = 'ls' and screen.index.method = 'ls'
# for stacking, screen.dat.method = 'all' and screen.index.method = 'ls'
res.screen <- screen.dat.index(Y = Y_vaccine, X = X_vaccine, X_markers = X_markers, fit.set = 'BAMA', 
                               screen.dat.method = 'all', screen.index.method = 'ls', obsWeights = weights_vaccine)
dat.X <- res.screen$dat
screen.var <- res.screen$screen.index.var



# 3. Fitting codes #
pred.mat <- matrix( NA, nrow = 10, ncol = 5 )
colnames( pred.mat ) <- can.set <- c('None','BAMA','Tcells','BAMA_Tcells','All')

fit.method <- 'ST'   # 'RF' for Random Forest analysis and 'ST' for Stacking analysis
fit.set <- 'BAMA'    # one of candidate sets: 'None', 'BAMA', 'Tcells', 'BAMA_Tcells', 'All'

for( i in 1:10 ){
  print(i)
  seed <- i
  
  if( fit.method == 'RF' ){
    # RF #
    res.temp <- get.rf.cvauc( dat = dat.X, cv.scheme = '5fold', obsWeights = weights_vaccine, method = 'sRF', seed = seed )
  } else if( fit.method == 'ST' ){
    # Stacking #
    res.temp <- get.st.cvauc( dat = dat.X, cv.scheme = '5fold', obsWeights = weights_vaccine, 
                              var.index = list(glm  = c(TRUE, rep(TRUE,3), rep(FALSE,ncol(dat.X$case)-3)), # baseline covariates
                                               rf   = c(TRUE, screen.var)),                                 # baseline covariates + screened variables
                                               #rf.1 = c(TRUE, rep(TRUE,ncol(dat.X$case)))),                # All variables
                              method = 'method.NNloglik', seed = seed )
  }
  
  pred.mat[i, fit.set] <- mean( res.temp )
}
apply(pred.mat, 2, mean)

write.table( pred.mat, file = '/Users/shan/Desktop/Paper/YFong/2.HVTN/Result/sRF_ulr_no.txt', col.names = T, sep = ',' )
