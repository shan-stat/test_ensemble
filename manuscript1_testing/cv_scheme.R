source( 'helper_cvauc.R' )
numCores <- as.numeric( Sys.getenv( 'SLURM_CPUS_ON_NODE' ) )



### alpha & beta ###
sep <- 0
n1 <- 25 ; n0 <- 125
#true.cvauc <- 0.5

### Generating result vector ###
res.mat <- matrix( NA, nrow = 1000, ncol = 4 )
rownames( res.mat ) <- 1:1000 ; colnames( res.mat ) <- c( 'LPO', '5fold', '10x5fold', 'random5fold' )

for( i in 1:1000 ){
  seed <- i
  dat <- sim.1( n1 = n1, n0 = n0, seed = seed, sep = sep )
  
  # LPO, 5-fold, Random5fold #
  res.mat[seed, '5fold'] <- get.cv.auc( dat = dat, cv.scheme = '5fold', method = 'MLR', numCores = numCores, seed = seed )
  
  # 10x5fold CV #
  #res.vec <- c()
  #for( j in 1:10 ){
  #  inner.seed <- j
  #  res.vec[j] <- get.cv.auc( dat = dat, cv.scheme = 'LPO', method = 'MLR', numCores = numCores, seed = inner.seed )
  #}
  #res.mat[seed, '10x5fold'] <- mean(res.vec)
}
res.mat

write.table( res.mat, file = paste0("n25_5fold.txt"), col.names = TRUE, sep = ',' )