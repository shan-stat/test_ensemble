source( 'helper_cvauc.R' )
numCores <- as.numeric( Sys.getenv( 'SLURM_CPUS_ON_NODE' ) )


### Setting method ###
mtd <- 'ULR'
alpha <- -0.7 ; beta <- 0.5 # MTCT
#alpha <- -1.9 ; beta <- c( 0.9, 0.7, 0.7, 0.5, -1 ) # RV144
beta.z.1 <- 0 ; beta.z.2 <- 0


### n1 & n0 ###
n1 <- 80 ; n0 <- 160 # MTCT
#n1 <- 40 ; n0 <- 200 # RV144
#n1 <- 120 ; n0 <- 120 # Balanced setting


### Generating data ###
seed <- as.numeric( Sys.getenv( 'SLURM_ARRAY_TASK_ID' ) )
dat <- sim.mtct( n1 = n1, n0 = n0, seed = seed, alpha = alpha, beta = beta, beta.z.1 = beta.z.1, beta.z.2 = beta.z.2 )
#dat <- sim.rv144( n1 = n1, n0 = n0, seed = seed, alpha = alpha, betas = beta, beta.z.1 = beta.z.1, beta.z.2 = beta.z.2 )


### Real n1 & n0 ###
n1 <- sum( dat$y ) ; n0 <- ( nrow( dat ) - n1 )


### Getting estimates ###
res.mat <- matrix( NA, nrow = 1000, ncol = 2 )
rownames( res.mat ) <- 1:1000 ; colnames( res.mat ) <- c( 'est.auc', mtd )

# Without covariates #
dat <- dat[, c('y', paste('x',1:(ncol(dat)-3),sep=''))] 
dat.temp <- list( case = subset(dat, y==1, select=-y), control = subset(dat, y==0, select=-y) )
dat.case.control <- rbind( dat.temp$case, dat.temp$control ) # For permutation distribution
res.mat[seed, 'est.auc'] <- get.cv.auc.ULR( dat = dat.temp, cv.scheme = 'nocv', numCores = numCores, seed = seed ) # nocv is for AUC 


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
    ref.distr[i] <- get.cv.auc.ULR( dat = dat.perm, cv.scheme = 'nocv', numCores = numCores, seed = i )
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
    ref.distr[i] <- get.cv.auc.ULR( dat = dat.perm, cv.scheme = 'nocv', numCores = numCores, seed = i )
  }
} else stop("wrong p.method")


### Calculating size & power ###
pvalue.test <- sum( res.mat[seed, 'est.auc'] < ref.distr ) / nperm
ifelse( pvalue.test < 0.05, res.mat[seed, mtd] <- 1, res.mat[seed, mtd] <- 0 )

write.table( res.mat[seed, ], file = paste0(mtd, "_out_", seed, ".txt"), col.names = TRUE, sep = ',' )