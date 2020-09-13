source( 'helper_cvauc.R' )


### Setting method ###
mtd <- 'BH'
alpha <- -0.7 ; beta <- 0.5 # MTCT
#alpha <- -1.9 ; beta <- c( 0.9, 0.7, 0.7, 0.5, -1 ) # RV144
beta.z.1 <- 0 ; beta.z.2 <- 0


### n1 & n0 ###
n1 <- 80 ; n0 <- 160 # MTCT
#n1 <- 40 ; n0 <- 200 # RV144
#n1 <- 120 ; n0 <- 120 # Balanced setting


### Getting estimates ###
res.mat <- matrix( NA, nrow = 1000, ncol = 2 )
rownames( res.mat ) <- 1:1000 ; colnames( res.mat ) <- c( 'est.pvalue', 'BH' )

for( i in 1:1000 ){
  seed <- i
  dat <- sim.mtct( n1 = n1, n0 = n0, seed = seed, alpha = alpha, beta = beta, beta.z.1 = beta.z.1, beta.z.2 = beta.z.2 )
  #dat <- sim.rv144( n1 = n1, n0 = n0, seed = seed, alpha = alpha, betas = beta, beta.z.1 = beta.z.1, beta.z.2 = beta.z.2 )
  
  # Calculating BH pvalue #
  pvalues <- c()
  
  for( j in 1:(ncol(dat)-3) ){
    dat.temp <- dat[, c('y','z1','z2',paste('x',j,sep=''))]
    set.seed(123)
    fit.temp <- glm( factor(y) ~ ., dat.temp, family=binomial ) 
    pvalues[j] <- summary(fit.temp)$coefficients[paste('x',j,sep=''), "Pr(>|z|)"]
  }
  res.temp <- { sort( pvalues ) < seq( (1/(ncol(dat)-3))*0.05, 0.05, length = (ncol(dat)-3) ) }
  
  if( sum(res.temp) == 0 ) {
    res.mat[i, 'est.pvalue'] <- 1 ; res.mat[i, 'BH'] <- 0
  } else {
    res.mat[i, 'est.pvalue'] <- pvalues[ max( which( res.temp ) ) ]
    ifelse( res.mat[i, 'est.pvalue'] < 0.05, res.mat[i, 'BH'] <- 1, res.mat[i, 'BH'] <- 0 )
  }
}

write.table( res.mat, file = paste( mtd, '_out.txt', sep = '' ), row.names = T, col.names = T, sep = ',' )