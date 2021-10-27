##### Section 2.1 #####
### Import helper function ###
source("helper_cvauc.R")

### 1. Setting arguments ###
Args <- commandArgs(trailingOnly=TRUE) 
if (length(Args)==0) {
  # batch.size and batch.number are arguments for running a batch script
  # sim.setting: "sim.1_n25_0" for simulation studies in Section 2 (you can change sample size from n25 to others)
  # basic.cv.setting: "5-fold" by default
  # proj: "cv_" for comparing cv schemes (e.g. "cv_LPO", "cv_5fold", "cv_10x5fold", and "cv_50xrandom4:1")
  Args=c(batch.size="1000",batch.number="10",sim.setting="sim.1_n25_0",proj="cv_50xrandom4:1")
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
  } else {
      
  }stop("wrong simple.size")
  
  # Simulated dataset #
  if (sim.model=="sim.1") {
    dat=sim.1( n1 = n1, n0 = n0, seed = seed, sep = sep )
  } else stop ("wrong sim.model")
  
  # 2-1. Comparison CV schemes #
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
  
  out
})
res
