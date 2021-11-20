##### Section 2.1 and Supplementary Materials Section C #####
### Import helper function ###
source("helper_cvauc.R")

### 1. Setting arguments ###
Args <- commandArgs(trailingOnly=TRUE) 
if (length(Args)==0) {
  # batch.size and batch.number are arguments for running a batch script
  # sim.setting: "sim.1_n25_0" for simulation studies in Section 2 (you can change sample size from n25 to others)
  # basic.cv.setting: "5-fold" by default
  # proj: "cv_" for comparing cv schemes (e.g. "cv_LPO", "cv_5fold", "cv_10x5fold", and "cv_50xrandom4:1") in Section 2.1
  #     : "perm" for permutation test; "LeDell" for LeDell CI-based test; "Benkeser" for Benkeser CI-based test in Supplementary Materials Section C
  Args=c(batch.size="1000",batch.number="10",sim.setting="sim.1_n25_0",proj="cv_5fold")
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
i=i+1; proj=Args[i]
verbose=ifelse(unix(),0,2)

nperm=1e3 # permutation replicates
cv.scheme='5fold' # basic cv scheme

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
  
  # Simulated dataset
  if (sim.model=="sim.1") {
    dat=sim.1( n1 = n1, n0 = n0, seed = seed, sep = sep )
  } else stop ("wrong sim.model")
  
  # Section 2.1 (Choice of CV schemes)
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
  
  # Supplementary Materials Section C (Size of confidence interval-based hypothesis testing approaches)
  if(proj=="perm"){
    # Permutation-based test
    est.cvauc <- get.cv.auc(dat=dat, cv.scheme=cv.scheme, seed=seed)
    # Permutation process
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
    
    # P-value
    pvalue.test <- sum(est.cvauc <= ref.distr) / nperm
    ifelse(pvalue.test < 0.05, reject <- 1, reject <- 0)
    out=c(est=est.cvauc, reject=reject)
    
  } else if(proj=="LeDell"){
    # LeDell CI-based test
    est.cvauc <- get.cv.auc(dat=dat, cv.scheme=cv.scheme, seed=seed)
    LeDell.ci <- get.cv.auc.LeDell(dat=dat, cv.scheme=cv.scheme, seed=1)
    ifelse( (0.5 < LeDell.ci['lb']) , reject <- 1, reject <- 0 ) # for type 1 error rate
    out=c(est=est.cvauc, reject=reject)
    
  } else if(proj=="Benkeser"){
    # Benkeser CI-based test
    est.cvauc <- get.cv.auc(dat=dat, cv.scheme=cv.scheme, seed=seed)
    cvtmle.ci <- get.cvtmle(dat=dat, seed=1)
    ifelse( (0.5 < cvtmle.ci['lb']) , reject <- 1, reject <- 0 ) # for type 1 error rate
    out=c(est=est.cvauc, reject=reject)
  }
  
  out
})
res
