source( 'helper_cvauc.R' )
numCores <- as.numeric( Sys.getenv( 'SLURM_CPUS_ON_NODE' ) )



### alpha & beta ###
sep <- 0
n1 <- 25 ; n0 <- 125
true.cvauc <- 0.5


### Generating result vector ###
res.mat <- matrix( NA, nrow = 1000, ncol = 3 )
rownames( res.mat ) <- 1:1000 ; colnames( res.mat ) <- c( 'est.cvauc', 'perc', 'basic' )
#res.mat <- matrix( NA, nrow = 1000, ncol = 5 )
#rownames( res.mat ) <- 1:1000 ; colnames( res.mat ) <- c( 'est.cvauc', 'perc', 'basic', 'perc_width', 'basic_width' )


seed <- as.numeric( Sys.getenv( 'SLURM_ARRAY_TASK_ID' ) )
dat <- sim.1( n1 = n1, n0 = n0, seed = seed, sep = sep )
res.mat[seed, 'est.cvauc'] <- get.cv.auc( dat = dat, cv.scheme = '5fold', method = 'MLR', numCores = numCores, seed = seed )


# Bootstrap CI #
res.cvauc <- rep( NA, length = 1000 )

for( k in 1:1000 ){
  print(k)
  # Generating bootstrapped sample #
  set.seed( k )
  idx.case <- sample( 1:nrow(dat$case), replace = TRUE )
  idx.control <- sample( 1:nrow(dat$control), replace = TRUE )
  boot.dat <- list( case = dat$case[idx.case, ], control = dat$control[idx.control, ] )
  
  # Calculating CV-AUC #
  res.cvauc[k] <- get.cv.auc( dat = boot.dat, cv.scheme = '5fold', method = 'MLR', numCores = numCores, seed = k )
}

# Percentile bootstrap CI #
perc.ci <- quantile( res.cvauc, c( 0.025, 0.975 ) )
perc.ci.lb <- perc.ci[1] ; perc.ci.ub <- perc.ci[2]
ifelse( (true.cvauc < perc.ci.lb ) , res.mat[seed, 'perc'] <- 1, res.mat[seed, 'perc'] <- 0 ) # Type 1 error
#ifelse( (perc.ci.lb <= true.cvauc) & (true.cvauc <= perc.ci.ub), res.mat[seed, 'perc'] <- 1, res.mat[seed, 'perc'] <- 0 ) # Coverage probability
#res.mat[seed, 'perc_width'] <- (perc.ci.ub - perc.ci.lb) # getting CI width

# Basic bootstrap CI #
basic.ci.lb <- 2*res.mat[seed, 'est.cvauc'] - quantile( res.cvauc, 0.975 )
basic.ci.ub <- 2*res.mat[seed, 'est.cvauc'] - quantile( res.cvauc, 0.025 )
ifelse( (true.cvauc < basic.ci.lb) , res.mat[seed, 'basic'] <- 1, res.mat[seed, 'basic'] <- 0 ) # Type 1 error
#ifelse( (basic.ci.lb <= true.cvauc) & (true.cvauc <= basic.ci.ub), res.mat[seed, 'basic'] <- 1, res.mat[seed, 'basic'] <- 0 ) # Coverage probability
#res.mat[seed, 'basic_width'] <- (basic.ci.ub - basic.ci.lb) # getting CI width

# Symmetric bootstrap CI #
#symm.ci.lb <- res.mat[seed, 'est.cvauc'] - quantile( abs( res.cvauc - res.mat[seed, 'est.cvauc'] ), 0.975 )
#symm.ci.ub <- res.mat[seed, 'est.cvauc'] + quantile( abs( res.cvauc - res.mat[seed, 'est.cvauc'] ), 0.975 )
#ifelse( (true.cvauc < symm.ci.lb) , res.mat[seed, 'symm'] <- 1, res.mat[seed, 'symm'] <- 0 )


write.table( res.mat[seed, ], file = paste0("boot_out_", seed, ".txt"), col.names = TRUE, sep = ',' )