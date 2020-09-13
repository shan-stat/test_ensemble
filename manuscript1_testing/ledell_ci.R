source( 'helper_cvauc.R' )
numCores <- as.numeric( Sys.getenv( 'SLURM_CPUS_ON_NODE' ) )



### Setting method ###
mtd <- 'LeDell'
sep <- 0
n1 <- 25 ; n0 <- 125
true.cvauc <- 0.5


### Generating result vector ###
res.mat <- matrix( NA, nrow = 1000, ncol = 2)
rownames( res.mat ) <- 1:1000 ; colnames( res.mat ) <- c( 'est.cvauc', mtd )
#res.mat <- matrix( NA, nrow = 1000, ncol = 3)
#rownames( res.mat ) <- 1:1000 ; colnames( res.mat ) <- c( 'est.cvauc', mtd, 'LeDell_width' )


for( i in 1:1000 ){
  seed <- i
  dat <- sim.1( n1 = n1, n0 = n0, seed = seed, sep = sep )
  res.mat[seed, 'est.cvauc'] <- get.cv.auc( dat = dat, cv.scheme = '5fold', method = 'MLR', numCores = numCores, seed = seed )
  
  # LeDell CI #
  LeDell.ci <- get.cv.auc.LeDell( dat = dat, cv.scheme = '5fold', seed = 1 )
  ifelse( (true.cvauc < LeDell.ci['lb']) , res.mat[seed, mtd] <- 1, res.mat[seed, mtd] <- 0 ) # Type 1 error
  #ifelse( (LeDell.ci['lb'] < true.cvauc) & (true.cvauc < LeDell.ci['ub']) , res.mat[seed, mtd] <- 1, res.mat[seed, mtd] <- 0 ) # Coverage probability
  #res.mat[seed, 'LeDell_width'] <- (LeDell.ci['ub'] - LeDell.ci['lb']) # getting CI width
}

write.table( res.mat, file = paste0(mtd, "_out.txt"), row.names = TRUE, col.names = TRUE, sep = ',' )