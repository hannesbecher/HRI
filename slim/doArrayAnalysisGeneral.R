setwd("~/git_repos/HRI/slim/")
source("analyseArraysUtils.R")
set.seed(123345)
#anaFun <- function(repID, l=2500){
anaFun <- function(repID, simID){
  if(missing(simID)) stop("Need to supply simulation ID.")
  
  parsLines <- readLines(paste0("params", simID, "/paramsV", simID))
  print(parsLines)
  eval(str2expression(parsLines))  
  
  gt2 <- read.gt(paste0(x = paste0("~/temp/HRI/v", simID, "/sim"), repID, ".gt2"))
  gt4 <- read.gt(paste0(x = paste0("~/temp/HRI/v", simID, "/sim"), repID, ".gt4"))
  
  
  samp2 <- lapply(c(10, 20, 40, 80), function(n) sampleGt(gt2, n))
  samp4 <- lapply(c(10, 20, 40, 80), function(n) sampleGt(gt4, n))
  
  # lists of SFSs
  sfs2 <- lapply(samp2, function(x) gt2sfs(x,"sfsS"))
  sfs4 <- lapply(samp4, function(x) gt2sfs(x,"sfsN"))

  # get conditional SFS proportions (deleterious)
  sfs2condRel <- lapply(sfs2, function(x){
    v <- x[2:(length(x)-1)]
    v / sum(v)
  })


  # 
  tp2 <- sapply(sfs2, theta_pi, persite=F)/LL
  names(tp2) <- c("tpS_10", "tpS_20", "tpS_40", "tpS_80")
  tw2 <- sapply(sfs2, theta_w, persite=F)/LL
  names(tw2) <- c("twS_10", "twS_20", "twS_40", "twS_80")
  ta2 <- sapply(sfs2, tajd)
  names(ta2) <- c("taS_10", "taS_20", "taS_40", "taS_80")
  de2 <- sapply(sfs2, deltaTheta)
  names(de2) <- c("deS_10", "deS_20", "deS_40", "deS_80")
  
  dtwp2 <- sapply(sfs2, deltaThetaPrime)
  names(dtwp2) <- c("dtwpS_10", "dtwpS_20", "dtwpS_40", "dtwpS_80")
  
  tp4 <- sapply(sfs4, theta_pi, persite=F)/LL
  names(tp4) <- c("tpN_10", "tpN_20", "tpN_40", "tpN_80")
  tw4 <- sapply(sfs4, theta_w, persite=F)/LL
  names(tw4) <- c("twN_10", "twN_20", "twN_40", "twN_80")
  ta4 <- sapply(sfs4, tajd)
  names(ta4) <- c("taN_10", "taN_20", "taN_40", "taN_80")
  de4 <- sapply(sfs4, deltaTheta)
  names(de4) <- c("deN_10", "deN_20", "deN_40", "deN_80")
  
  dtwp4 <- sapply(sfs4, deltaThetaPrime)
  names(dtwp4) <- c("dtwpN_10", "dtwpN_20", "dtwpN_40", "dtwpN_80")
  
  pBar2 <- 1-getPBar(gt2, l=LL)
  pBar4 <- 1-getPBar(gt4, l=LL)
  
  pBar2No0 <- 1-sum(sfs2[[4]])/LL # analog to older SLiM version
  
  gt2v <- gt2[apply(gt2, 1, sd) !=0,]
  gt4v <- gt4[apply(gt4, 1, sd) !=0,]
  n2 <- nrow(gt2v)
  n4 <- nrow(gt4v)
  
  pwLD2 <- pairwiseLD(gt2)
  LDmean2 <- mean(pwLD2)
  pqMean2 <- mean(pairwisePQProd(gt2))
  minus2LdS <- -2*LDmean2*n2
  
  
  pwLD4 <- pairwiseLD(gt4)
  LDmean4 <- mean(pwLD4)
  pqMean4 <- mean(pairwisePQProd(gt4))
  minus2LdN <- -2*LDmean4*n4
  
  
  gHat <- log(pBar2/(1-pBar2)) # selected sites, scalar
  gammaExp <- -2*NN*ss
  B <- gHat/gammaExp
  Bprime <- unname(tp4[4])/gammaExp# neutral sites, computed from pi, only for SFS80
  
  
  dataRow <- c(LL=LL, ss=ss, NN=NN, uu=uu, kk=kk, gammaExp=gammaExp,
               tp2, tw2, ta2, de2,dtwp2,
               tp4, tw4, ta4, de4,dtwp4,
               pBar2=pBar2, pBar4=pBar4,
               LD2=LDmean2/sqrt(pqMean2), minus2LdS=minus2LdS,
               LD4=LDmean4/sqrt(pqMean4), minus2LdN=minus2LdN,
               n2=n2, n4=n4,
               B=B, Bprime=Bprime,
               do.call(c, sfs2),
               do.call(c, sfs4)
               )
  
  return(dataRow)
}
r0 <- anaFun("0001", "02")
#r0



library(parallel)
t0 <- Sys.time()
a <- mclapply(1:10, function(x) anaFun(sprintf("%04d", x), "02"), mc.cores = 10)
Sys.time() - t0
# about 4 mins on Gemmaling






# Per-run rows --------------------------------------------------------------------------------
b <- as.data.frame(do.call(rbind, a))
b



meanSE(b[,1])
names(b)
# mSe <- t(apply(b[,1:46], 2, meanSE))
# 
# mSe[42,]
# mSe[-c(36, 38),]
# 
#t(apply(b, 2, meanSE))[1:32,]
# pBar <- 0.672
# pBar <- mSe[33,1] # check position!
# gammaHat <- log(pBar/(1-pBar))
# gamma <- 20 # depends in sim parameters!
# gammaHat/gamma
# getwd()
# write.table(b,
#             "allNumsv03_240610.tsv",
#             row.names = F,
#             quote = F)
# write.table(mSe,
#             "mSe_v03_240610.tsv",
#             col.names = NA,
#             row.names = T,
#             quote = F)
# 


# Avg SFS -------------------------------------------------------------------------------------

sfsPart <- b[,startsWith(names(b), "sfs")]

# get averaged conditional SFS
propCondSfs <- function(x, fro, to){
  ll <- to-fro+1
  y <- x[,fro:to]
  z <- colSums(y)
  z[1] <- 0
  z[ll] <- 0
  sfs <- z/sum(z)
  sfs[2:(ll-1)]
}
head(sfsPart)
cSfsS10 <- propCondSfs(sfsPart, 1,11)
cSfsS20 <- propCondSfs(sfsPart, 12,32)
cSfsS40 <- propCondSfs(sfsPart, 33,73)
cSfsS80 <- propCondSfs(sfsPart, 74,154)

cSfsN10 <- propCondSfs(sfsPart, 155,165)
cSfsN20 <- propCondSfs(sfsPart, 166,186)
cSfsN40 <- propCondSfs(sfsPart, 187,227)
cSfsN80 <- propCondSfs(sfsPart, 228,308)




# Stats from PK -------------------------------------------------------------------------------


l10 <- readLines("v03sfs10.txt")
l20 <- readLines("v03sfs20.txt")
l40 <- readLines("v03sfs40.txt")
l80 <- readLines("v03sfs80.txt")



E10 <- sapply(l10, secCol, USE.NAMES = F)
E20 <- sapply(l20, secCol, USE.NAMES = F)
E40 <- sapply(l40, secCol, USE.NAMES = F)
E80 <- sapply(l80, secCol, USE.NAMES = F)

# P&K expectations as list
sfsPK <- list(E10, E20, E40, E80)


pic <- theta_pi(sfsPK[[1]])
pic
theta_w(sfsPK[[4]])


# pi = theta_w * a_n * pic
# p_seg = theta_w * a_n () # i.e. depends on sample size

# delta_theta_w = 1 - pi/theta_w
# = 1 - theta_w*a_n*pic/theta_w
# = 1 - a_n*pic

# delta-theta-w = 1 – 1/(pic-c x an)
1 - theta_pi(c(0, sfsPK[[1]], 0))*a_sub_1(10)
1 - theta_pi(c(0, sfsPK[[2]], 0))*a_sub_1(20)
1 - theta_pi(c(0, sfsPK[[3]], 0))*a_sub_1(40)
1 - theta_pi(c(0, sfsPK[[4]], 0))*a_sub_1(80)




# delta theta w prime -------------------------------------------------------------------------
# 
# sapply(c(10, 20, 40, 80), deltaThetaMax)
# deltaTheta(c(0, cSfsN10, 0))
# deltaTheta(c(0, cSfsN20, 0))
# deltaTheta(c(0, cSfsN40, 0))
# deltaTheta(c(0, cSfsN80, 0))
deltaThetaPrime(c(0, cSfsN10, 0))
deltaThetaPrime(c(0, cSfsN20, 0))
deltaThetaPrime(c(0, cSfsN40, 0))
deltaThetaPrime(c(0, cSfsN80, 0))


deltaThetaPrime(c(0, sfsPK[[1]], 0))
deltaThetaPrime(c(0, sfsPK[[2]], 0))
deltaThetaPrime(c(0, sfsPK[[3]], 0))
deltaThetaPrime(c(0, sfsPK[[4]], 0))

# get conditional SFS proportions (neutral sites)
sfs4condRel <- lapply(sfs4, function(x){
  v <- x[2:(length(x)-1)]
  v / sum(v)
})

# Delta theta_w prime
dtwp <- sapply(sfs4condRel, deltaThetaPrime)
names(dtwp) <- c("dtwp10", "dtwp20", "dtwp40", "dtwp80")



# Stuff ---------------------------------------------------------------------------------------
sss0 <- b[1,274:354]
sss1 <- sss0; sss1[1] <- 0; sss1[81] <- 0
sss2 <- sss1/sum(sss1)
tajd(sss0)
tajd(sss1)
tajd(sss2)



curve(dchisq(x,2),-1,10)
curve(pchisq(x,2),-1,10)


# Wrap it all up ####
# analyses of reps
# get CIs
# extract avg SFS (N/S)
# 
library(parallel)
t0 <- Sys.time()

Sys.time() - t0

simID <- c("02")

analyseAll <- function(simID, maxRep=10, nCores=10){
  # TODO autodetermine number of runs
  a <- mclapply(1:maxRep, function(x) anaFun(sprintf("%04d", x), "02"), mc.cores = 10)
  b <- as.data.frame(do.call(rbind, a))
  statCI <- t(apply(b, 2, meanSE))
  
  # SFS 
  
  sfsPart <- b[,startsWith(names(b), "sfs")]
  
  # get averaged conditional SFS
  propCondSfs <- function(x, fro, to){
    ll <- to-fro+1
    y <- x[,fro:to]
    z <- colSums(y)
    z[1] <- 0
    z[ll] <- 0
    sfs <- z/sum(z)
    sfs[2:(ll-1)]
  }
  head(sfsPart)
  cSfsS10 <- propCondSfs(sfsPart, 1,11)
  cSfsS20 <- propCondSfs(sfsPart, 12,32)
  cSfsS40 <- propCondSfs(sfsPart, 33,73)
  cSfsS80 <- propCondSfs(sfsPart, 74,154)
  
  cSfsN10 <- propCondSfs(sfsPart, 155,165)
  cSfsN20 <- propCondSfs(sfsPart, 166,186)
  cSfsN40 <- propCondSfs(sfsPart, 187,227)
  cSfsN80 <- propCondSfs(sfsPart, 228,308)
  
  
  l10 <- readLines(paste0("params", simID, "/v", simID, "sfs10.txt"))
  l20 <- readLines(paste0("params", simID, "/v", simID, "sfs20.txt"))
  l40 <- readLines(paste0("params", simID, "/v", simID, "sfs40.txt"))
  l80 <- readLines(paste0("params", simID, "/v", simID, "sfs80.txt"))
  
  
  
  E10 <- sapply(l10, secCol, USE.NAMES = F)
  E20 <- sapply(l20, secCol, USE.NAMES = F)
  E40 <- sapply(l40, secCol, USE.NAMES = F)
  E80 <- sapply(l80, secCol, USE.NAMES = F)
  
  # P&K expectations as list
  # need to add padding 0s because P&K give only segregating probs
  sfsPK <- list(c(0, E10, 0), c(0, E20, 0), c(0, E40, 0), c(0, E80, 0)) 
  
  
  pic <- theta_pi(sfsPK[[1]])
  pic
  theta_w_from_avg <- theta_w(sfsPK[[4]])
  

  
  dtwp10PK <- deltaThetaPrime(sfsPK[[1]])
  dtwp20PK <- deltaThetaPrime(sfsPK[[2]])
  dtwp40PK <- deltaThetaPrime(sfsPK[[3]])
  dtwp80PK <- deltaThetaPrime(sfsPK[[4]])
  
  list(mSE = statCI,
       
       dtwp10PK=dtwp10PK,
       dtwp20PK=dtwp20PK,
       dtwp40PK=dtwp40PK,
       dtwp80PK=dtwp80PK,
       
       sfsPart=sfsPart # check if this works correctly and whether required at all. SFS is mean SE object as well.
       )
  
}
aaa <- analyseAll("02", maxRep = 100)


propCondSfs(aaa$mSE[57:67,1])
propCondSfs(aaa$mSE, 57, 67)
