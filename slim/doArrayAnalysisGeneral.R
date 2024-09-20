# Use this script to compute stats for sets of simulations

setwd("~/git_repos/HRI/slim/")
source("analyseArraysUtils.R")
#set.seed(123345)
#anaFun <- function(repID, l=2500){
anaFun <- function(repID, simID){
  if(missing(simID)) stop("Need to supply simulation ID.")
  
  parsLines <- readLines(paste0("params", simID, "/paramsV", simID))
  print(parsLines)
  eval(str2expression(parsLines))  
  
  # genotype matrices (whole population) with position as "rownames"
  gt2 <- read.gt(paste0("params", simID, "/data/sim", repID, ".gt2"))
  gt4 <- read.gt(paste0("params", simID, "/data/sim", repID, ".gt4"))
  
  
  # lists of samples of 10/20/40/80 haploids
  samp2 <- lapply(c(10, 20, 40, 80), function(n) sampleGt(gt2, n))
  samp4 <- lapply(c(10, 20, 40, 80), function(n) sampleGt(gt4, n))
  
  # lists of SFSs corresponding to the samples taken above
  sfs2 <- lapply(samp2, function(x) gt2sfs(x,"sfsS"))
  sfs4 <- lapply(samp4, function(x) gt2sfs(x,"sfsN"))

  # get conditional SFS proportions (deleterious sites)
  sfs2condRel <- lapply(sfs2, function(x){
    v <- x[2:(length(x)-1)]
    v / sum(v)
  })
  
  # get individuals' fitnesesses
  indW <- exp(colSums(log(gt2 * ss + 1))) # individuals are in columns
  indWsum <- colSums(gt2 * ss) + 1 # individuals are in columns
  varW <- var(indW)
  varWsum<- var(indWsum)
  
  
  # pi
  tp2 <- sapply(sfs2, theta_pi, persite=F)/LL
  names(tp2) <- c("tpS_10", "tpS_20", "tpS_40", "tpS_80")
  
  # Watterson's th
  tw2 <- sapply(sfs2, theta_w, persite=F)/LL
  names(tw2) <- c("twS_10", "twS_20", "twS_40", "twS_80")
  
  # TajD
  ta2 <- sapply(sfs2, tajd)
  names(ta2) <- c("taS_10", "taS_20", "taS_40", "taS_80")
  
  # Delta theta
  de2 <- sapply(sfs2, deltaTheta)
  names(de2) <- c("deS_10", "deS_20", "deS_40", "deS_80")
  
  # Delta theta prime (= del the / del the max)
  dtwp2 <- sapply(sfs2, deltaThetaPrime)
  names(dtwp2) <- c("dtwpS_10", "dtwpS_20", "dtwpS_40", "dtwpS_80")
  
  # same for neutral sites
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
  dtwm4 <- sapply(sfs4, function(x) deltaThetaMax(length(x)-1))
  names(dtwm4) <- c("dtwmN_10", "dtwmN_20", "dtwmN_40", "dtwmN_80")
  
  # avg allele frequencies (deleterious site val is used to estimate B)
  # the Li-Bulmer eqtn give p bar for the beneficial allele, so need to to 1-freq here:
  pBar2 <- 1-getPBar(gt2, l=LL)
  pBar4 <- 1-getPBar(gt4, l=LL)
  
  # what's this, again? legacy?
  #pBar2No0 <- 1-sum(sfs2[[4]])/LL # analog to older SLiM version
  
  # variant position only
  gt2v <- gt2[apply(gt2, 1, sd) !=0,]
  gt4v <- gt4[apply(gt4, 1, sd) !=0,]
  # How many variant sites in each class? (del and neutral)
  n2 <- nrow(gt2v)
  n4 <- nrow(gt4v)
  
  # pi vals for variant deleterious sites
  sitePis2 <- getSitePis(gt2v)
  varWexp <- sum(sitePis2)/2*ss^2
  
  
  pwLD2 <- pairwiseLD(gt2) # a vector (lower triangle of pairwise comparisons)
  LDmean2 <- mean(pwLD2)
  LDsumSel <- sum(pwLD2)
  pqMean2 <- mean(pairwisePQProd(gt2))
  minus2LdS <- -2*LDmean2*n2
  
  
  pwLD4 <- pairwiseLD(gt4)
  LDmean4 <- mean(pwLD4)  # adjust to all sites! bit only variant ones needed for correlation based
  LDsumNeu <- sum(pwLD4)
  pqMean4 <- mean(pairwisePQProd(gt4))
  minus2LdN <- -2*LDmean4*n4
  
  
  gHat <- log(pBar2/(1-pBar2)) # selected sites, scalar
  gammaExp <- -2*NN*ss
  B <- gHat/gammaExp
  Bprime <- unname(tp4[4])/2/NN/uu# neutral sites, computed from observed pi and expected w/o sle, only for SFS80
  
  
  dataRow <- c(LL=LL, ss=ss, NN=NN, uu=uu, kk=kk, gammaExp=gammaExp,
               tp2, tw2, ta2, de2,dtwp2,
               tp4, tw4, ta4, de4,dtwp4,
               pBar2=pBar2, pBar4=pBar4,
               #
               LDmeanSelVar=LDmean2, LDsumSel=LDsumSel, LDmean2ss=LDmean2 * abs(ss), LD2=LDmean2/sqrt(pqMean2),# minus2LdS=minus2LdS,
               LDmeanNeuVar=LDmean4, LDsumNeu=LDsumNeu, LDmean4ss=LDmean4 * abs(ss), LD4=LDmean4/sqrt(pqMean4),# minus2LdN=minus2LdN,
               n2=n2, n4=n4,
               B=B, Bprime=Bprime,
               varW=varW, varWexp=varWexp, varWsum=varWsum,
               do.call(c, sfs2),
               do.call(c, sfs4)
               )
  
  return(dataRow)
}




r0 <- anaFun("0002", "06")
#r0[1:60]

# 
# 
# library(parallel)
# t0 <- Sys.time()
# a <- mclapply(1:10, function(x) anaFun(sprintf("%04d", x), "02"), mc.cores = 10)
# Sys.time() - t0
# about 4 mins on Gemmaling






# Per-run rows --------------------------------------------------------------------------------
# b <- as.data.frame(do.call(rbind, a))
# b



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
# 
# sfsPart <- b[,startsWith(names(b), "sfs")]
# 
# # get averaged conditional SFS
# propCondSfs <- function(x, fro, to){
#   ll <- to-fro+1
#   y <- x[,fro:to]
#   z <- colSums(y)
#   z[1] <- 0
#   z[ll] <- 0
#   sfs <- z/sum(z)
#   sfs[2:(ll-1)]
# }
# head(sfsPart)
# cSfsS10 <- propCondSfs(sfsPart, 1,11)
# cSfsS20 <- propCondSfs(sfsPart, 12,32)
# cSfsS40 <- propCondSfs(sfsPart, 33,73)
# cSfsS80 <- propCondSfs(sfsPart, 74,154)
# 
# cSfsN10 <- propCondSfs(sfsPart, 155,165)
# cSfsN20 <- propCondSfs(sfsPart, 166,186)
# cSfsN40 <- propCondSfs(sfsPart, 187,227)
# cSfsN80 <- propCondSfs(sfsPart, 228,308)
# 
# 
# 
# 
# # Stats from PK -------------------------------------------------------------------------------
# 
# 
# l10 <- readLines("v03sfs10.txt")
# l20 <- readLines("v03sfs20.txt")
# l40 <- readLines("v03sfs40.txt")
# l80 <- readLines("v03sfs80.txt")
# 
# 
# 
# E10 <- sapply(l10, secCol, USE.NAMES = F)
# E20 <- sapply(l20, secCol, USE.NAMES = F)
# E40 <- sapply(l40, secCol, USE.NAMES = F)
# E80 <- sapply(l80, secCol, USE.NAMES = F)
# 
# # P&K expectations as list
# sfsPK <- list(E10, E20, E40, E80)
# 
# 
# pic <- theta_pi(sfsPK[[1]])
# pic
# theta_w(sfsPK[[4]])
# 
# 
# # pi = theta_w * a_n * pic
# # p_seg = theta_w * a_n () # i.e. depends on sample size
# 
# # delta_theta_w = 1 - pi/theta_w
# # = 1 - theta_w*a_n*pic/theta_w
# # = 1 - a_n*pic
# 
# # delta-theta-w = 1 – 1/(pic-c x an)
# 1 - theta_pi(c(0, sfsPK[[1]], 0))*a_sub_1(10)
# 1 - theta_pi(c(0, sfsPK[[2]], 0))*a_sub_1(20)
# 1 - theta_pi(c(0, sfsPK[[3]], 0))*a_sub_1(40)
# 1 - theta_pi(c(0, sfsPK[[4]], 0))*a_sub_1(80)
# 
# 
# 
# 
# # delta theta w prime -------------------------------------------------------------------------
# # 
# # sapply(c(10, 20, 40, 80), deltaThetaMax)
# # deltaTheta(c(0, cSfsN10, 0))
# # deltaTheta(c(0, cSfsN20, 0))
# # deltaTheta(c(0, cSfsN40, 0))
# # deltaTheta(c(0, cSfsN80, 0))
# deltaThetaPrime(c(0, cSfsN10, 0))
# deltaThetaPrime(c(0, cSfsN20, 0))
# deltaThetaPrime(c(0, cSfsN40, 0))
# deltaThetaPrime(c(0, cSfsN80, 0))
# 
# 
# deltaThetaPrime(c(0, sfsPK[[1]], 0))
# deltaThetaPrime(c(0, sfsPK[[2]], 0))
# deltaThetaPrime(c(0, sfsPK[[3]], 0))
# deltaThetaPrime(c(0, sfsPK[[4]], 0))
# 
# # get conditional SFS proportions (neutral sites)
# sfs4condRel <- lapply(sfs4, function(x){
#   v <- x[2:(length(x)-1)]
#   v / sum(v)
# })
# 
# # Delta theta_w prime
# dtwp <- sapply(sfs4condRel, deltaThetaPrime)
# names(dtwp) <- c("dtwp10", "dtwp20", "dtwp40", "dtwp80")
# 
# 
# 
# # Stuff ---------------------------------------------------------------------------------------
# sss0 <- b[1,274:354]
# sss1 <- sss0; sss1[1] <- 0; sss1[81] <- 0
# sss2 <- sss1/sum(sss1)
# tajd(sss0)
# tajd(sss1)
# tajd(sss2)
# 
# 
# 
# curve(dchisq(x,2),-1,10)
# curve(pchisq(x,2),-1,10)
# 

# Wrap it all up ####
# per parameter set
## analyses of reps
### get CIs on these
### extract avg SFS (N/S)
## P&K SFSs
## table Delta theta_w prime for P&K (no CI) and from sim with CI (neutral sites)
## table col LD for qp and -2, both with CI
## table col pBar, B, and B' from sim and analytic/P&K


library(parallel)
#simID <- c("02")

analyseAll <- function(simID, maxRep=10, nCores=10){
  # TODO autodetermine number of runs
  a <- mclapply(1:maxRep, function(x) anaFun(sprintf("%04d", x), simID), mc.cores = nCores)
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
  #head(sfsPart)
  cSfsS10 <- propCondSfs(sfsPart, 1,11)
  cSfsS20 <- propCondSfs(sfsPart, 12,32)
  cSfsS40 <- propCondSfs(sfsPart, 33,73)
  cSfsS80 <- propCondSfs(sfsPart, 74,154)
  
  cSfsN10 <- propCondSfs(sfsPart, 155,165)
  cSfsN20 <- propCondSfs(sfsPart, 166,186)
  cSfsN40 <- propCondSfs(sfsPart, 187,227)
  cSfsN80 <- propCondSfs(sfsPart, 228,308)
  
  # Get P&K SFS data
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
       ldNums = statCI[rownames(statCI) %in% c("LD2", "LDmeanSelVar", "LDsumSel", "LD4", "LDmeanNeuVar", "LDsumNeu"),],
       pbbNums = statCI[rownames(statCI) %in% c("pBar2", "B", "Bprime"),],
       dtwp = statCI[rownames(statCI) %in% c("dtwpN_10", "dtwpN_20", "dtwpN_40", "dtwpN_80"),],
       varNums = statCI[rownames(statCI) %in% c("varW", "varWexp", "varWsum"),],
       sfsNums = statCI[startsWith(rownames(statCI), "sfs"),],
       sfsPK=sfsPK,
       dtwp10PK=dtwp10PK,
       dtwp20PK=dtwp20PK,
       dtwp40PK=dtwp40PK,
       dtwp80PK=dtwp80PK,
       cSfsN = list(c(0, cSfsN10, 0), c(0, cSfsN20, 0), c(0, cSfsN40, 0), c(0, cSfsN80, 0)),
       cSfsS = list(c(0, cSfsS10, 0), c(0, cSfsS20, 0), c(0, cSfsS40, 0), c(0, cSfsS80, 0))
       #sfsPart=sfsPart # check if this works correctly and whether required at all. SFS is mean SE object as well.
       )
  
}
a02 <- analyseAll("02", maxRep = 200, nCores = 10) 
#a02 <- analyseAll("02", maxRep = 10, nCores = 10)
a02$ldNums
a02$varNums
a03 <- analyseAll("03", maxRep = 200, nCores = 10) # there are 200 replicates
a04 <- analyseAll("04", maxRep = 200, nCores = 10) # L=10k, takes long!, there are 200 replicates
a05 <- analyseAll("05", maxRep = 200, nCores = 10) # L=10k, takes long!, there are 200 replicates
# need to add P&K sfs files before this can run:
t0 <- Sys.time()
a06 <- analyseAll("06", maxRep = 200, nCores = 10)
Sys.time() - t0
rownames(aaa$mSE)
#rownames(bbb)
a06[[1]]

# make code to combine mean and CIs into one col!
array2tableCol <- function(arr, sigfig=3){
  rn <- rownames(arr)
  vals <- apply(signif(arr, sigfig), 1, function(x) paste0(x[1], " (CI: ", x[2], ", ", x[3], ")"))
  #  names(vals) <- rn
  vals
}


array2tableCol(aaa$ldNums, 4)
array2tableCol(aaa$pbbNums, 4)
array2tableCol(aaa$dtwp, 4)



propCondSfs(aaa$mSE[57:67,1])
propCondSfs(aaa$mSE, 57, 67)



# Combine param sets --------------------------------------------------------------------------
a02$mSE[,1]
a02$sfsNums
a02$sfsPK

#ps <- list(a02, a03) # 10 each for debug
ps <- list(a02, a03, a04, a05)
#saveRDS(ps, "ps.rds")
ps <- readRDS("ps.rds")
length(ps[[1]])
length(a06)
ps[[1]][[1]]
ps6 <- c(ps, list(a06))
#ps <- ps6
#saveRDS(ps6, "ps6.rds")
#ps <- readRDS("ps.rds")


makeTable <- function(parSets, itemName){
  n <- length(parSets)
  nums <- lapply(parSets, function(x) x[[itemName]])
  #print(nums)
  tab <- sapply(nums, array2tableCol, sigfig=4)
  colnames(tab) <- paste("params", 2:(ncol(tab)+1), sep = "")
  tab
}
#apply(signif(sapply(ps, function(x) c(x$dtwp10PK, x$dtwp20PK, x$dtwp40PK, x$dtwp80PK)), 4), 2, as.character)
makeDtwpTable <- function(parSets){
  fromSim <- makeTable(parSets, "dtwp")
  fromPK <- apply(signif(sapply(parSets, function(x) c(x$dtwp10PK, x$dtwp20PK, x$dtwp40PK, x$dtwp80PK)), 4), 2, as.character)
  colnames(fromPK) <- paste("params", 2:(length(parSets)+1), "PK", sep = "")
  #print(fromSim)
  #print(fromPK)
  cols <- lapply(1:length(parSets), function(x) cbind(fromSim[,x], fromPK[,x]))
  tab <- do.call(cbind,  cols)
  colnames(tab) <- as.vector(rbind(colnames(fromSim), colnames(fromPK)))
  return(tab)
}


getSfsN <- function(parSets){
  msfs <- lapply(parSets, function(x) do.call(c, x$cSfsN))
  do.call(cbind, msfs)
}

getPkSfs <- function(parSets){
  sfs <- lapply(parSets, function(x) do.call(c, x$sfsPK))
  do.call(cbind, sfs)
}

makeSfsTab <- function(parSets){
  sfsN <- getSfsN(ps)
  pkSfsN <- getPkSfs(ps)
  cols <- lapply(1:length(parSets), function(x) cbind(sfsN[,x], pkSfsN[,x]))
  tab <- do.call(cbind,  cols)
  colnames(tab) <- paste0(c("sfsSim", "sfsPK"), rep(2:(length(parSets)+1), each=2))
  signif(tab, 4)
}


pN <- t(makeTable(ps, "pbbNums"))
lN <- t(makeTable(ps, "ldNums"))
vN <- t(makeTable(ps, "varNums"))
sN <- makeSfsTab(ps)
dN <- t(makeDtwpTable(ps))
gN <- makeTable(ps, "mSE")[1:46,] # adjust range




library(openxlsx2)

wb <- wb_workbook(creator = "HB") %>%
  wb_add_worksheet(sheet = "general") %>%
  wb_add_data(x = data.frame(gN), row_names = TRUE) %>%
  wb_add_worksheet(sheet = "fitness variance") %>%
  wb_add_data(x = data.frame(vN), row_names = TRUE) %>%
  wb_add_worksheet(sheet = "LD") %>%
  wb_add_data(x = data.frame(lN), row_names = TRUE) %>%
  wb_add_worksheet(sheet = "B") %>%
  wb_add_data(x = data.frame(pN), row_names = TRUE) %>%
  wb_add_worksheet(sheet = "Dtwp") %>%
  wb_add_data(x = data.frame(dN), row_names = TRUE) %>%
  wb_add_worksheet(sheet = "SFS (neutral)") %>%
  wb_add_data(x = data.frame(sN), row_names = TRUE)


for(i in wb_get_sheet_names(wb)){
  wb_set_col_widths(wb, sheet = i, cols=1:100, widths="auto")  
}

#wb_save(wb, "Simulations-12.xlsx")


