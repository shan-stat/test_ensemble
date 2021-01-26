### Import helper function ###
setwd("set your local directory where helper_cvauc.R is located")
source("helper_cvauc.R")



### 1. Setting process arguments ###
Args <- commandArgs(trailingOnly=TRUE) 
if (length(Args)==0) {
  # batch.size and batch.number are need to a batch script
  # sim.setting: sim.1_n25_0 for simulation studies in Section 2 (you can change sample size of n25 listed in Section 2)
  # basic.cv.setting: 5-fold cross validation
  # proj: boot for bootstrap CI-based test, perm for permutation-based test, 
  #       Ledell for Ledell CI-based test, cv_ for cv scheme (you can choose one of cv_LPO, cv_5fold, cv_10x5fold, and cv_50xrandom4:1)
  Args=c(batch.size="1000",batch.number="1",sim.setting="sim.1_n80_0",basic.cv.setting="5fold",proj="boot")
}
myprint(Args)
i=0;
i=i+1; batch.size=as.numeric(Args[i])
i=i+1; batch=as.numeric(Args[i])
seeds=1:batch.size+batch.size*(batch-1); names(seeds)=seeds
myprint(batch, batch.size)
i=i+1; sim.setting=Args[i]; tmp=strsplit(sim.setting,"_")[[1]]
sim.model=tmp[1]
sample.size=tmp[2]
sep=as.numeric(tmp[3]) 
i=i+1; basic.cv.setting=Args[i]
cv.scheme=basic.cv.setting
i=i+1; proj=Args[i]

nperm=1e3 # permutation replicates
nboot=1e3 # bootstrap replicates
verbose=ifelse(unix(),0,2)
make.plot=FALSE # TRUE for figure 1
if(make.plot) {seeds=1:9; proj="plot"; mypdf(mfrow=c(3,3), file='boot_perm_dist')}



### 2. Experiments ###
res=sapply(seeds, simplify="array", function (seed) {
  
  myprint(seed)
  
  # Sample size #
  if (sample.size=="n25") {
    n1=25; n0=125
  } else if (sample.size=="n40") {
    n1=40; n0=200
  } else if (sample.size=="n80") {
    n1=80; n0=160
  } else if (sample.size=="n500") {
    n1=500; n0=500
  } else if (sample.size=="n2000") {
    n1=2000; n0=2000
  } else stop("wrong simple.size")
  
  # Simulated dataset #
  if (sim.model=="sim.1") {
    dat=sim.1( n1 = n1, n0 = n0, seed = seed, sep = sep )
  } else stop ("wrong sim.model")
  
  # 2-1. Comparison of CV schemes #
  if(startsWith(proj,"cv")){
    if(proj != 'cv_10x5fold'){
      # LPO, 5fold, 50xRandom4:1 #
      cv.mtd <- substr(proj,start=4,stop=nchar(proj))
      est.cvauc <- get.cv.auc(dat=dat, cv.scheme=cv.mtd, seed=seed)
    } else if(proj == 'cv_10x5fold'){
      # 10x5fold #
      res.vec <- c()
      for( j in 1:10 ){
        inner.seed <- j
        res.vec[j] <- get.cv.auc(dat=dat, cv.scheme="5fold", seed=inner.seed)
      }
      est.cvauc <- mean(res.vec)
    }
    out=c(est=est.cvauc)
  } 
  
  # 2-2. Fitting tests #
  if(proj=="perm"){
    # Permutation-based test #
    est.cvauc <- get.cv.auc(dat=dat, cv.scheme='5fold', seed=seed)
    
    # Permutation process #
    N <- n0 + n1
    dat.case.control <- rbind( dat$case, dat$control )
    ref.distr <- c()
    for( i in 1:nperm ){
      set.seed( i )
      ind.1 <- sample( 1:N, n1 )    # cases in the permutated dataset
      ind.0 <- setdiff( 1:N, ind.1 ) # controls in the permutated dataset
      dat.perm <- list()
      dat.perm$case <- dat.case.control[ind.1,]
      dat.perm$control <- dat.case.control[ind.0,]
      ref.distr[i] <- get.cv.auc(dat=dat.perm, cv.scheme=cv.scheme, seed=i)
    }
    
    # P-value #
    pvalue.test <- sum(est.cvauc < ref.distr) / nperm
    ifelse(pvalue.test < 0.05, reject <- 1, reject <- 0)
    out=c(est=est.cvauc, reject=reject)
    
  } else if(proj=="boot"){
    # Bootstrap CI-based test #
    est.cvauc <- get.cv.auc(dat=dat, cv.scheme=cv.scheme, seed=seed)
    
    # Bootstrap process #
    ref.distr <- c()
    for( k in 1:nboot ){
      set.seed( k )
      idx.case <- sample( 1:nrow(dat$case), replace = TRUE )
      idx.control <- sample( 1:nrow(dat$control), replace = TRUE )
      boot.dat <- list( case = dat$case[idx.case, ], control = dat$control[idx.control, ] )
      ref.distr[k] <- get.cv.auc(dat=boot.dat, cv.scheme=cv.scheme, seed=k)
    }
    
    # Percentile bootstrap CI #
    perc.ci <- quantile(ref.distr, c(0.05, 0.95)) # alpha=0.05; (1-2alpha)% CI
    perc.ci.lb <- perc.ci[1] ; perc.ci.ub <- perc.ci[2]
    ifelse( (0.5 < perc.ci.lb ) , reject <- 1, reject <- 0 )
    out=c(est=est.cvauc, reject=reject)

  } else if(proj=="Ledell"){
    # LeDell's CI-based test #
    est.cvauc <- get.cv.auc(dat=dat, cv.scheme=cv.scheme, seed=seed)
    
    # LeDell CI #
    LeDell.ci <- get.cv.auc.LeDell(dat=dat, cv.scheme=cv.scheme, seed=1)
    ifelse( (0.5 < LeDell.ci['lb']) , reject <- 1, reject <- 0 ) # 0.5 is thoretical cvauc under the null; Type 1 error
    out=c(est=est.cvauc, reject=reject)
  } 
  
  # 2-3. Figure #
  if(make.plot){
    # Estimated CV-AUC #
    est.cvauc <- get.cv.auc(dat=dat, cv.scheme='5fold', seed=seed)
    
    # Permutation distribution #
    N <- n0 + n1
    dat.case.control <- rbind( dat$case, dat$control )
    perm.distr <- c()
    for( i in 1:nperm ){
      set.seed( i )
      ind.1 <- sample( 1:N, n1 )    # cases in the permutated dataset
      ind.0 <- setdiff( 1:N, ind.1 ) # controls in the permutated dataset
      dat.perm <- list()
      dat.perm$case <- dat.case.control[ind.1,]
      dat.perm$control <- dat.case.control[ind.0,]
      perm.distr[i] <- get.cv.auc(dat=dat.perm, cv.scheme=cv.scheme, seed=i)
    }
    
    # Bootstrap process #
    boot.distr <- c()
    for( k in 1:nboot ){
      set.seed( k )
      idx.case <- sample( 1:nrow(dat$case), replace = TRUE )
      idx.control <- sample( 1:nrow(dat$control), replace = TRUE )
      boot.dat <- list( case = dat$case[idx.case, ], control = dat$control[idx.control, ] )
      boot.distr[k] <- get.cv.auc(dat=boot.dat, cv.scheme=cv.scheme, seed=k)
    }
  
    plot(density(perm.distr), main=seed, col='blue', xlim=c(0.2,0.8), ylim=c(0,8)); abline(v=est.cvauc, lty=2); abline(v=0.5)
    lines(density(boot.distr), col='red')
    out=c(est=est.cvauc)
  }

  out
})
if(make.plot) dev.off()
res
