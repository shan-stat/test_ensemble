source( 'helper_cvauc.R' )
numCores <- as.numeric( Sys.getenv( 'SLURM_CPUS_ON_NODE' ) )



### alpha & beta ###
sep <- 0
n1 <- 25 ; n0 <- 125
true.cvauc <- 0.5

### Generating result vector ###
res.mat <- matrix( NA, nrow = 1000, ncol = 2 )
rownames( res.mat ) <- 1:1000 ; colnames( res.mat ) <- c( 'est.cvauc', 'perm' )


seed <- as.numeric( Sys.getenv( 'SLURM_ARRAY_TASK_ID' ) )
dat <- sim.1( n1 = n1, n0 = n0, seed = seed, sep = sep )
res.mat[seed, 'est.cvauc'] <- get.cv.auc( dat = dat, cv.scheme = '5fold', method = 'MLR', numCores = numCores, seed = seed )
dat.case.control <- rbind( dat$case, dat$control ) # For fitting permutation distribution


### Permutation process ###
nperm <- 1000

N <- n0 + n1
if( choose(N, n1) <= nperm ) p.method <- "exact" else p.method <- "Monte Carlo"

if( p.method == "exact" ){
  ref.distr <- c()
  indexes <- combn( 1:N, n1 )
  for( i in 1:nperm ){
    ind.1 <- indexes[,i]
    ind.0 <- setdiff(1:N, ind.1)
    dat.perm <- list()
    dat.perm$case <- data.frame( dat.case.control[ind.1,] )
    dat.perm$control <- data.frame( dat.case.control[ind.0,] )
    ref.distr[i] <- get.cv.auc( dat = dat.perm, cv.scheme = '5fold', method = 'MLR', numCores = numCores, seed = i )
  }
} else if ( p.method == "Monte Carlo" ){
  ref.distr <- c()
  for( i in 1:nperm ){
    set.seed( i )
    
    # Without covariates #
    ind.1 <- sample( 1:N, n1 )    # cases in the permutated dataset
    ind.0 <- setdiff( 1:N, ind.1 ) # controls in the permutated dataset
    dat.perm <- list()
    dat.perm$case <- dat.case.control[ind.1,]
    dat.perm$control <- dat.case.control[ind.0,]
    ref.distr[i] <- get.cv.auc( dat = dat.perm, cv.scheme = '5fold', method = 'MLR', numCores = numCores, seed = i )
  }
} else stop("wrong p.method")


### Calculating size & power ###
pvalue.test <- sum( res.mat[seed, 'est.cvauc'] < ref.distr ) / nperm
ifelse( pvalue.test < 0.05, res.mat[seed, 'perm'] <- 1, res.mat[seed, 'perm'] <- 0 )

write.table( res.mat[seed, ], file = paste0("perm_out_", seed, ".txt"), col.names = TRUE, sep = ',' )