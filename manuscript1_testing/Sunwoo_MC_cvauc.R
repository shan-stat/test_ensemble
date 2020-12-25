source( 'helper_cvauc.R' )

# process arguments
Args <- commandArgs(trailingOnly=TRUE) 
if (length(Args)==0) {
  # proj: boot for bootstrap, perm for permutation test, Ledell for Ledell CI, cv for cross-validation scheme
  Args=c(batch.size="10",batch.number="1",sim.setting="sim.1_n25_0",basic.cv.setting="5fold",proj="cv_10x5fold")
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
seed=1 # temp seed for debugging use, does not matter

begin=Sys.time()
res=sapply(seeds, simplify="array", function (seed) {
  
  myprint(seed)
  t.0=Sys.time()   
  
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
  } else if (sample.size=="n10000") {
    n1=NULL; n0=NULL; n=1e4
  } else stop("wrong simple.size")
  
  if (sim.model=="sim.1") {
    dat=sim.1( n1 = n1, n0 = n0, seed = seed, sep = sep )
  } else stop ("wrong sim.model")
  
  # perform testing
  if(proj=="boot"){
    # Estimating CV-AUC #
    est.cvauc <- get.cv.auc(dat=dat, cv.scheme=cv.scheme, seed=seed)
    
    # Bootstrap CI #
    res.cvauc <- rep( NA, length = nboot )
    for( k in 1:nboot ){
      print(k)
      # Generating bootstrapped sample #
      set.seed( k )
      idx.case <- sample( 1:nrow(dat$case), replace = TRUE )
      idx.control <- sample( 1:nrow(dat$control), replace = TRUE )
      boot.dat <- list( case = dat$case[idx.case, ], control = dat$control[idx.control, ] )
      # Calculating CV-AUC #
      res.cvauc[k] <- get.cv.auc(dat=boot.dat, cv.scheme=cv.scheme, seed=k)
    }
    
    # Percentile bootstrap CI #
    perc.ci <- quantile(res.cvauc, c(0.025, 0.975))
    perc.ci.lb <- perc.ci[1] ; perc.ci.ub <- perc.ci[2]
    ifelse( (0.5 < perc.ci.lb ) , reject <- 1, reject <- 0 ) # 0.5 is thoretical cvauc under the null; Type 1 error
    out=c(est=est.cvauc, reject=reject)
    
  } else if(proj=="Ledell"){
    # Estimating CV-AUC #
    est.cvauc <- get.cv.auc(dat=dat, cv.scheme=cv.scheme, seed=seed)
    
    # LeDell CI #
    LeDell.ci <- get.cv.auc.LeDell(dat=dat, cv.scheme=cv.scheme, seed=1)
    ifelse( (0.5 < LeDell.ci['lb']) , reject <- 1, reject <- 0 ) # 0.5 is thoretical cvauc under the null; Type 1 error
    out=c(est=est.cvauc, reject=reject)
  } else if(proj=="perm"){
    # Estimating CV-AUC #
    est.cvauc <- get.cv.auc(dat=dat, cv.scheme='5fold', seed=seed)
    
    ### Permutation process ###
    N <- n0 + n1
    if( choose(N, n1) <= nperm ) p.method <- "exact" else p.method <- "Monte Carlo"
    dat.case.control <- rbind( dat$case, dat$control ) # For fitting permutation distribution
    if( p.method == "exact" ){
      ref.distr <- c()
      indexes <- combn( 1:N, n1 )
      for( i in 1:nperm ){
        ind.1 <- indexes[,i]
        ind.0 <- setdiff(1:N, ind.1)
        dat.perm <- list()
        dat.perm$case <- data.frame( dat.case.control[ind.1,] )
        dat.perm$control <- data.frame( dat.case.control[ind.0,] )
        ref.distr[i] <- get.cv.auc(dat=dat.perm, cv.scheme=cv.scheme, seed=i)
      }
    } else if ( p.method == "Monte Carlo" ){
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
    } else stop("wrong p.method")
    
    ### Permutation test ###
    pvalue.test <- sum(est.cvauc < ref.distr) / nperm
    ifelse(pvalue.test < 0.05, reject <- 1, reject <- 0)
    out=c(est=est.cvauc, reject=reject)
  } else if(startsWith(proj,"cv")){
    if(proj != 'cv_10x5fold'){
      # LPO, 5-fold, Random5fold #
      cv.mtd <- substr(proj,start=4,stop=nchar(proj))
      est.cvauc <- get.cv.auc(dat=dat, cv.scheme=cv.mtd, seed=seed)
    } else if(proj == 'cv_10x5fold'){
      # 10x5fold CV #
      res.vec <- c()
      for( j in 1:10 ){
        inner.seed <- j
        res.vec[j] <- get.cv.auc(dat=dat, cv.scheme="5fold", seed=inner.seed)
      }
      est.cvauc <- mean(res.vec)
    }
    out=c(est=est.cvauc)
  }
})
print(date()%.%". Time used: "%.%format(Sys.time()-begin))

# save res
foldername="res_"%.%proj%.%"/"; if(!file.exists(foldername)) dir.create(foldername)
foldername=foldername%.%sim.setting%.%"/"; if(!file.exists(foldername)) dir.create(foldername)
save(res, file=foldername%.%"/batch"%.%formatInt(batch, 3)%.%".Rdata")

