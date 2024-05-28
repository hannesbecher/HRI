setwd("~/git_repos/HRI/slim/params02/")
source("analyseArrays.R")
set.seed(123345)
anaFun <- function(repID, l=2500){
  gt2 <- read.gt(paste0(x = "~/temp/HRI/sim", repID, ".gt2"))
  gt4 <- read.gt(paste0(x = "~/temp/HRI/sim", repID, ".gt4"))
  
  
  samp2 <- lapply(c(10, 20, 40, 80), function(n) sampleGt(gt2, n))
  samp4 <- lapply(c(10, 20, 40, 80), function(n) sampleGt(gt4, n))
  
  sfs2 <- lapply(samp2, gt2sfs)
  sfs4 <- lapply(samp4, gt2sfs)

  sfs2condRel <- lapply(sfs2, function(x){
    v <- x[2:(length(x)-1)]
    v / sum(v)
  })

  sfs4condRel <- lapply(sfs4, function(x){
    v <- x[2:(length(x)-1)]
    v / sum(v)
  })

  cSq <- mapply(function(o,e){
    sum((o-e)^2/e)
  },
  sfs4condRel,
  sfsExp
  )
  # #print(cSq) # debug
  # 
  tp2 <- sapply(sfs2, theta_pi, persite=F)/l
  tw2 <- sapply(sfs2, theta_w, persite=F)/l
  ta2 <- sapply(sfs2, tajd)
  de2 <- sapply(sfs2, deltaTheta)
  
  
  
  tp4 <- sapply(sfs4, theta_pi, persite=F)/l
  tw4 <- sapply(sfs4, theta_w, persite=F)/l
  ta4 <- sapply(sfs4, tajd)
  de4 <- sapply(sfs4, deltaTheta)
  
  pBar2 <- 1-getPBar(gt2, l=l)
  pBar4 <- 1-getPBar(gt4, l=l)
  
  pBar2No0 <- 1-sum(sfs2[[4]])/l # analog to older SLiM version
  
  gt2v <- gt2[apply(gt2, 1, sd) !=0,]
  gt4v <- gt4[apply(gt4, 1, sd) !=0,]
  n2 <- nrow(gt2v)
  n4 <- nrow(gt4v)
  
  LDmean2 <- sum(pairwiseLD(gt2))
  rmean2 <- meanN(pairwiser(gt2), l*(l-1)/2)
  LDmean4 <- sum(pairwiseLD(gt4))
  rmean4 <- meanN(pairwiser(gt4), l*(l-1)/2)
  
  B <- log(pBar2/(1-pBar2)) # selected sites, scalar
  Bprime <- tp4/(2*1000*1e-5)# neutral sites, computed from pi, vector
  
  dataRow <- c(tp2, tw2, ta2, de2,
               tp4, tw4, ta4, de4,
               pBar2, pBar4,
               LDmean2, rmean2, LDmean4, rmean4,
               n2, n4,
               B, Bprime,
               cSq
               )
  names(dataRow) <- c("tpS_10", "tpS_20", "tpS_40", "tpS_80",
                      "twS_10", "twS_20", "twS_40", "twS_80",
                      "taS_10", "taS_20", "taS_40", "taS_80",
                      "deS_10", "deS_20", "deS_40", "deS_80",
                      
                      "tpN_10", "tpN_20", "tpN_40", "tpN_80",
                      "twN_10", "twN_20", "twN_40", "twN_80",
                      "taN_10", "taN_20", "taN_40", "taN_80",
                      "deN_10", "deN_20", "deN_40", "deN_80",
                      
                      "pBar2", "pBar4",
                      
                      "ldSumS", "rmeanS", "ldSumN", "rmeanN",
                      "n2", "n4",
                      "B", "bPrime10", "bPrime20", "bPrime40", "bPrime80",
                      "cSq10", "cSq20", "cSq40", "cSq80"
                      )
  dataRow
}
#r0 <- anaFun("0001")
#r0



library(parallel)
t0 <- Sys.time()
a <- mclapply(1:100, function(x) anaFun(sprintf("%04d", x)), mc.cores = 10)
Sys.time() - t0
b <- as.data.frame(do.call(rbind, a))



meanSE <- function(vec){
  m <- mean(vec)
  se <- sd(vec) / sqrt(length(vec))
  c(mean=m, lo = m + qnorm(0.025)*se, hi = m + qnorm(0.975)*se)
}

meanSE(b[,1])
mSe <- t(apply(b, 2, meanSE))

mSe[-c(36, 38),]

#t(apply(b, 2, meanSE))[1:32,]
log(7.092162e-01/(1-7.092162e-01))/2
# getwd()
# write.table(b,
#             "allNums.tsv",
#             row.names = F,
#             quote = F)
# write.table(mSe,
#             "mSe.tsv",
#             col.names = NA,
#             row.names = T,
#             quote = F)


# LD things ---------------------------------------------------------------

plot(b$ldSumS / b$n2 /(b$n2-1) * 2 * -2 * b$n2,

b$ldSumS/ 1250/2499 * -5000
)
grid()

ldSvec <- b$ldSumS / b$n2 /(b$n2-1) * 2 * -2 * b$n2
ldNvec <- b$ldSumN / b$n4 /(b$n4-1) * 2 * -2 * b$n4
meanSE(ldSvec)
meanSE(ldNvec)


# Chi squared things ------------------------------------------------------




sfsExp
a_sub_1 = function(n) {
  sum(1/(1:(n-1)))
}
sum(1/(1:9*a_sub_1(10)))

expNeut <- function(n){
  1/(1:(n-1)*a_sub_1(n))
}
cSq <- function(e, o){
  sum((e-o)^2/e)
}
cSq(expNeut(10), sfsExp[[1]])
cSq(expNeut(20), sfsExp[[2]])
cSq(expNeut(40), sfsExp[[3]])
cSq(expNeut(80), sfsExp[[4]])


# from sims

sfsFun <- function(repID, l=2500){
  gt4 <- read.gt(paste0(x = "~/temp/HRI/sim", repID, ".gt4"))
  
  samp4 <- lapply(c(10, 20, 40, 80), function(n) sampleGt(gt4, n))
  sfs4 <- lapply(samp4, gt2sfs)
  sfs4var <- lapply(sfs4, function(x) x[2:(length(x)-1)])
  sfs4varRel <- lapply(sfs4var, function(x) x/sum(x))
  do.call(c, sfs4varRel)
}
sfsFun("0002")
ss <- mclapply(1:100, function(x) sfsFun(sprintf("%04d", x)), mc.cores = 10)
tt <- do.call(rbind, ss)
uu <- colMeans(tt)
plot(uu)
vv <- list(uu[1:9],
           uu[10:28],
           uu[29:67],
           uu[68:146])
mapply(function(o,e){
  sum((o-e)^2/e)
},
vv,
sfsExp
)
sfsPK
sfsExp[[1]]
