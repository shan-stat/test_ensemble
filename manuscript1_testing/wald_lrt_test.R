source( 'helper_cvauc.R' )


### alpha & beta ###
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
res.mat1 <- matrix( NA, nrow = 1000, ncol = 2 )
rownames( res.mat1 ) <- 1:1000 ; colnames( res.mat1 ) <- c( 'est.pvalue', 'Wald' )
res.mat2 <- matrix( NA, nrow = 1000, ncol = 2 )
rownames( res.mat2 ) <- 1:1000 ; colnames( res.mat2 ) <- c( 'est.pvalue', 'LRT' )


# Without covariates
dat <- dat[, c('y', paste('x',1:(ncol(dat)-3),sep=''))] 
dat.temp <- list( case = subset(dat, y==1, select=-y), control = subset(dat, y==0, select=-y) )
dat.case.control <- rbind( dat.temp$case, dat.temp$control ) # For permutation distribution
res.temp <- get.pvalue.test( dat = dat )
res.mat1[seed, 'est.pvalue'] <- res.temp[1] # Wald test 
res.mat2[seed, 'est.pvalue'] <- res.temp[2] # LRT

### Permutation process ###
nperm <- 1000

N <- n0 + n1
if( choose(N, n1) <= nperm ) p.method <- "exact" else p.method <- "Monte Carlo"

if( p.method == "exact" ){
  ref.distr1 <- c()
  ref.distr2 <- c()
  indexes <- combn( 1:N, n1 )
  for( i in 1:nperm ){
    ind.1 <- indexes[,i]
    ind.0 <- setdiff(1:N, ind.1)
    dat.perm <- list()
    dat.perm$case <- data.frame( y=1, dat.case.control[ind.1,] )
    dat.perm$control <- data.frame( y=0, dat.case.control[ind.0,] )
    dat.perm <- rbind( dat.perm$case, dat.perm$control )
    res.perm <- get.pvalue.test( dat = dat.perm )
    ref.distr1[i] <- res.perm[1] # Wald test 
    ref.distr2[i] <- res.perm[2] # LRT
  }
} else if ( p.method == "Monte Carlo" ){
  ref.distr1 <- c()
  ref.distr2 <- c()
  for( i in 1:nperm ){
    set.seed( i )
    
    # Without covariates #
    ind.1 <- sample( 1:N, n1 )    # cases in the permutated dataset
    ind.0 <- setdiff( 1:N, ind.1 ) # controls in the permutated dataset
    dat.perm <- list()
    dat.perm$case <- data.frame( y=1, dat.case.control[ind.1,] )
    dat.perm$control <- data.frame( y=0, dat.case.control[ind.0,] )
    dat.perm <- rbind( dat.perm$case, dat.perm$control )
    
    res.perm <- get.pvalue.test( dat = dat.perm )
    ref.distr1[i] <- res.perm[1] # Wald test 
    ref.distr2[i] <- res.perm[2] # LRT
  }
} else stop("wrong p.method")


### Calculating size & power ###
pvalue.test1 <- sum( ref.distr1 < res.mat1[seed, 'est.pvalue'] ) / nperm
ifelse( pvalue.test1 < 0.05, res.mat1[seed, 'Wald'] <- 1, res.mat1[seed, 'Wald'] <- 0 )
pvalue.test2 <- sum( ref.distr2 < res.mat2[seed, 'est.pvalue'] ) / nperm
ifelse( pvalue.test2 < 0.05, res.mat2[seed, 'LRT'] <- 1, res.mat2[seed, 'LRT'] <- 0 )

write.table( res.mat1[seed, ], file = paste0("Wald_out_", seed, ".txt"), col.names = TRUE, sep = ',' )
write.table( res.mat2[seed, ], file = paste0("LRT_out_", seed, ".txt"), col.names = TRUE, sep = ',' )