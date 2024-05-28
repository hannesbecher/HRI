setwd("~/git_repos/HRI/slim/params02/")
source("analyseArrays.R")
set.seed(123345)
anaFun <- function(repID, l=2500){
  gt2 <- read.gt(paste0(x = "~/temp/HRI/sim", repID, ".gt2"))
  gt4 <- read.gt(paste0(x = "~/temp/HRI/sim", repID, ".gt4"))
  
  
  samp2 <- lapply(c(10, 20, 40, 80), function(n) sampleGt(gt2, n))
  samp4 <- lapply(c(10, 20, 40, 80), function(n) sampleGt(gt4, n))
  
  # lists of SFSs
  sfs2 <- lapply(samp2, function(x) gt2sfs(x,"S"))
  sfs4 <- lapply(samp4, function(x) gt2sfs(x,"N"))

  # get conditional SFS proportions (deleterious)
  sfs2condRel <- lapply(sfs2, function(x){
    v <- x[2:(length(x)-1)]
    v / sum(v)
  })

  # get conditional SFS proportions (neutral sites)
  sfs4condRel <- lapply(sfs4, function(x){
    v <- x[2:(length(x)-1)]
    v / sum(v)
  })

  # Brian's "chi squared" per replicate
  cSq <- mapply(function(o,e){
    sum((o-e)^2/e)
  },
  sfs4condRel,
  sfsExp
  )
  names(cSq) <- c("cSq10", "cSq20", "cSq40", "cSq80")
  # #print(cSq) # debug
  
  # 
  tp2 <- sapply(sfs2, theta_pi, persite=F)/l
  names(tp2) <- c("tpS_10", "tpS_20", "tpS_40", "tpS_80")
  tw2 <- sapply(sfs2, theta_w, persite=F)/l
  names(tw2) <- c("twS_10", "twS_20", "twS_40", "twS_80")
  ta2 <- sapply(sfs2, tajd)
  names(ta2) <- c("taS_10", "taS_20", "taS_40", "taS_80")
  de2 <- sapply(sfs2, deltaTheta)
  names(de2) <- c("deS_10", "deS_20", "deS_40", "deS_80")
  
  
  
  tp4 <- sapply(sfs4, theta_pi, persite=F)/l
  names(tp4) <- c("tpN_10", "tpN_20", "tpN_40", "tpN_80")
  tw4 <- sapply(sfs4, theta_w, persite=F)/l
  names(tw4) <- c("twN_10", "twN_20", "twN_40", "twN_80")
  ta4 <- sapply(sfs4, tajd)
  names(ta4) <- c("taN_10", "taN_20", "taN_40", "taN_80")
  de4 <- sapply(sfs4, deltaTheta)
  names(de4) <- c("deN_10", "deN_20", "deN_40", "deN_80")
  
  pBar2 <- 1-getPBar(gt2, l=l)
  pBar4 <- 1-getPBar(gt4, l=l)
  
  pBar2No0 <- 1-sum(sfs2[[4]])/l # analog to older SLiM version
  
  gt2v <- gt2[apply(gt2, 1, sd) !=0,]
  gt4v <- gt4[apply(gt4, 1, sd) !=0,]
  n2 <- nrow(gt2v)
  n4 <- nrow(gt4v)
  
  LDmean2 <- mean(pairwiseLD(gt2))
  pqMean2 <- mean(pairwisePQProd(gt2))
  rmean2 <- meanN(pairwiser(gt2), l*(l-1)/2)
  LDmean4 <- mean(pairwiseLD(gt4))
  pqMean4 <- mean(pairwisePQProd(gt4))
  rmean4 <- meanN(pairwiser(gt4), l*(l-1)/2)
  
  B <- log(pBar2/(1-pBar2)) # selected sites, scalar
  Bprime <- unname(tp4[4])/(2*1000*1e-5)# neutral sites, computed from pi, only for SFS80
  
  
  dataRow <- c(tp2, tw2, ta2, de2,
               tp4, tw4, ta4, de4,
               pBar2=pBar2, pBar4=pBar4,
               LD2=LDmean2/sqrt(pqMean2), rmean2=rmean2, LD4=LDmean4/sqrt(pqMean4), rmean4=rmean4,
               n2=n2, n4=n4,
               B=B, Bprime=Bprime,
               cSq,
               do.call(c, sfs2),
               do.call(c, sfs4)
               )
  # names(dataRow) <- c(
  #                    
  #                     
  #                     "ldSumS", "rmeanS", "ldSumN", "rmeanN",
  #                     "n2", "n4",
  #                     "B", "bPrime10", "bPrime20", "bPrime40", "bPrime80",
  #                     "cSq10", "cSq20", "cSq40", "cSq80"
  #                     )
  dataRow
}
#r0 <- anaFun("0001")
#r0



library(parallel)
t0 <- Sys.time()
a <- mclapply(1:100, function(x) anaFun(sprintf("%04d", x)), mc.cores = 10)
Sys.time() - t0
# about 4 mins on Gemmaling


# Per-run rows --------------------------------------------------------------------------------
b <- as.data.frame(do.call(rbind, a))
b


meanSE <- function(vec){
  m <- mean(vec)
  se <- sd(vec) / sqrt(length(vec))
  c(mean=m, CIlo = m + qnorm(0.025)*se, CIhi = m + qnorm(0.975)*se)
}

meanSE(b[,1])
names(b)
mSe <- t(apply(b[,1:46], 2, meanSE))

mSe[42,]
mSe[-c(36, 38),]

#t(apply(b, 2, meanSE))[1:32,]
# pBar <- 0.797
pBar <- mSe[33,1] # check position!
gammaHat <- log(pBar/(1-pBar))
gamma <- 2 # depends in sim parameters!
gammaHat/gamma
# getwd()
write.table(b,
            "allNums240528.tsv",
            row.names = F,
            quote = F)
write.table(mSe,
            "mSe240528.tsv",
            col.names = NA,
            row.names = T,
            quote = F)


# Avg SFS -------------------------------------------------------------------------------------

sfsPart <- b[,47:354]

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
# same as computed in analyseArrays.R
cSq(expNeut(10), sfsExp[[1]])
cSq(expNeut(20), sfsExp[[2]])
cSq(expNeut(40), sfsExp[[3]])
cSq(expNeut(80), sfsExp[[4]])

cSq(expNeut(80), sfsPK[[4]])

mapply(function(x,y) cSq(x, y),
       sfsPK, 
       list(expNeut(10), expNeut(20), expNeut(40), expNeut(80))
)
# simuations (neutral sites)
mapply(function(x,y) cSq(x, y),
       list(cSfsN10,cSfsN20,cSfsN40,cSfsN80), 
       list(expNeut(10), expNeut(20), expNeut(40), expNeut(80))
)

barplot(rbind(cSfsN10, expNeut(10)), beside=T)
barplot(rbind(cSfsN20, expNeut(20)), beside=T)
barplot(rbind(cSfsN80, expNeut(80)), beside=T)









# Stats from PK -------------------------------------------------------------------------------
tajd(sfsPK[[4]])
pic <- theta_pi(c(0, sfsPK[[1]], 0))
pic
theta_w(c(0, sfsPK[[4]], 0))
tajd(c(0, sfsPK[[4]], 0))

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




# Stuff ---------------------------------------------------------------------------------------
sss0 <- b[1,274:354]
sss1 <- sss0; sss1[1] <- 0; sss1[81] <- 0
sss2 <- sss1/sum(sss1)
tajd(sss0)
tajd(sss1)
tajd(sss2)
