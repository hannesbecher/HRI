
source("analyseArrays.R")

anaFun <- function(repID, l=2500){
  gt2 <- read.gt(paste0(x = "~/temp/HRI/sim", repID, ".gt2"))
  gt4 <- read.gt(paste0(x = "~/temp/HRI/sim", repID, ".gt4"))
  
  
  samp2 <- lapply(c(10, 20, 40, 80), function(n) sampleGt(gt2, n))
  samp4 <- lapply(c(10, 20, 40, 80), function(n) sampleGt(gt4, n))
  
  sfs2 <- lapply(samp2, gt2sfs)
  sfs4 <- lapply(samp4, gt2sfs)
  
  tp2 <- sapply(sfs2, theta_pi, persite=F)/l
  tw2 <- sapply(sfs2, theta_w, persite=F)/l
  ta2 <- sapply(sfs2, tajd)
  de2 <- sapply(sfs2, deltaTheta)
  
  tp4 <- sapply(sfs4, theta_pi, persite=F)/l
  tw4 <- sapply(sfs4, theta_w, persite=F)/l
  ta4 <- sapply(sfs4, tajd)
  de4 <- sapply(sfs4, deltaTheta)
  
  
  LDmean2 <- meanN(pairwiseLD(gt2), l*(l-1)/2)
  rmean2 <- meanN(pairwiser(gt2), l*(l-1)/2)
  LDmean4 <- meanN(pairwiseLD(gt4), l*(l-1)/2)
  rmean4 <- meanN(pairwiser(gt4), l*(l-1)/2)
  #plot(pairwiseLD(gt2), pairwiser2(gt2))
  
  
  dataRow <- c(tp2, tw2, ta2, de2,
               tp4, tw4, ta4, de4,
               LDmean2, rmean2, LDmean4, rmean4)
  names(dataRow) <- c("tpS_10", "tpS_20", "tpS_40", "tpS_80",
                      "twS_10", "twS_20", "twS_40", "twS_80",
                      "taS_10", "taS_20", "taS_40", "taS_80",
                      "deS_10", "deS_20", "deS_40", "deS_80",
                      
                      "tpN_10", "tpN_20", "tpN_40", "tpN_80",
                      "twN_10", "twN_20", "twN_40", "twN_80",
                      "taN_10", "taN_20", "taN_40", "taN_80",
                      "deN_10", "deN_20", "deN_40", "deN_80",
                      
                      "LDmeanS", "rmeanS", "LDmeanN", "rmeanN")
  dataRow
}
r0 <- anaFun("0001")
r0
sprintf("%03d", 0)
library(parallel)

a <- mclapply(1:100, function(x) anaFun(sprintf("%04d", x)), mc.cores = 10)
b <- as.data.frame(do.call(rbind, a))

boxplot(b[,1:4])
names(b)
boxplot(b[,33:36])
meanSE <- function(vec){
  m <- mean(vec)
  se <- sd(vec) / sqrt(length(vec))
  c(mean=m, lo = m - qnorm(0.025)*se, hi = m + qnorm(0.025)*se)
}
meanSE(b[,1])
t(apply(b, 2, meanSE))
t(apply(b, 2, meanSE))[1:32,]
