# this is a script of utility functions, to be used bu doArrayAnalysisGeneral.R


# import gt files
# assumes 1st col is position, all other are haploid GTs (0/1)
read.gt <- function(x){
  a <-read.table(x, sep=",")
  a <- a[order(a[,1]),]
  gtM <- as.matrix(a[,2:ncol(a)])
  dimnames(gtM) <- list(a[,1],
                     paste0("i", sprintf("%03d",0:(ncol(a)-2)))
                     )
  
  gtM
}



gt2sfs <- function(gts,pref){
  # gts is assumed to be a 2D 0/1 array with the row names containing positions
  n <- ncol(gts)
  cts <- rowSums(gts)
  sfs <- sapply(0:n, function(k){
    sum(cts==k)
  })
  names(sfs) <- paste0(pref,n,"_",0:n)
  sfs
}

# GTs are sampled from haplotype matrix (by sampling columns)
sampleGt <- function(gt, n){
  # returns a 2D array
  gt[, sample(x = 1:ncol(gt), size = n, replace = F)]
}

# Stats functions --------------------------------------------------------------
theta_pi = function(sfs, persite=F) {
  n = length(sfs) - 1
  if(n <= 1) { return(NA) }
  x = n:0/n
  tot = sum(sfs)
  if(tot==0) { return(NA) }
  if(persite) {
    (2*n)/(n-1) * sum(sfs*x*(1 - x)) / tot
  } else {
    (2*n)/(n-1) * sum(sfs*x*(1 - x))
  }
}

theta_w = function(sfs, persite=F) {
  n = length(sfs) - 1
  if(n <= 1) { return(NA) }
  tot = sum(sfs)
  if(tot==0) { return(NA) }
  if(persite) {
    seg_sites(sfs)/(a_sub_1(n) * tot)
  } else {
    seg_sites(sfs)/a_sub_1(n)
  }
}

# 
# return the average allele frequency at a site.
# note this is called "theta" by bcftools
# 
mean_allele_freq = function(sfs) {
  n = length(sfs) - 1
  sum(0:n * sfs)/(n*sum(sfs))
}

# 
# return denominator of Tajima's D statistic for given n and s
# n = no. alleles sampled
# s = no. segregating sites
# 
tajd_denom = function(n, s) {
  a_1 = a_sub_1(n)
  a_2 = a_sub_2(n)
  b_1 = (n+1)/(3*(n-1))
  b_2 = (2*(n^2 + n + 3)) / (9*n*(n-1))
  e_1 = (b_1 - 1/a_1) / a_1
  e_2 = (b_2 - (n+2)/(a_1*n) + a_2/(a_1^2)) / (a_1^2 + a_2) 
  sqrt(e_1*s + e_2*s*(s-1))
}


# 
# estimate Tajima's D statistic for given sfs
# note that the denominator is 0 when no. segregating sites is 0
# so the statistic is undefined!
# 
seg_sites = function(sfs) {
  sum(sfs[2:(length(sfs)-1)])
}
# 
# a_sub_1 =  1 + 1/2 + 1/3 + ... + 1/(n-1)
# 
a_sub_1 = function(n) {
  sum(1/(1:(n-1)))
}


# 
# a_sub_2 =  1 + 1/(2^2) + 1/(3^2) + ... + 1/((n-1)^2)
# 
a_sub_2 = function(n) {
  sum(1/((1:(n-1))^2))
}
tajd = function(sfs) {
  n = length(sfs) - 1
  s = seg_sites(sfs)
  (theta_pi(sfs) - theta_w(sfs)) / tajd_denom(n, s)
}

deltaTheta <- function(sfs, persite=F){
  1 - theta_pi(sfs, persite=persite)/theta_w(sfs, persite=persite)
}
deltaThetaMax <- function(n){
  1 - 2*a_sub_1(n)/n
}

deltaThetaPrime <- function(sfs, persite=T){
  deltaTheta(sfs, persite = persite) / deltaThetaMax(length(sfs)-1)
}


# LD for a 2-row matrix
ldD <- function(gtRows){
  # takes a 2-row haplotype array
  
  l <- ncol(gtRows) # number of haploid inds
  
  # haplotype frequencies
  xx <- proportions(sapply(0:3, function(x) sum(x == c(1,2) %*% gtRows)))
  
  
  # allele frequencies at the loci
  qA <- sum(gtRows[1,])/l
  qB <- sum(gtRows[2,])/l
  
  
  xx[4] - qA*qB
}

pqProd <- function(gtRows){
  # takes a 2-row haplotype array
  
  l <- ncol(gtRows) # number of haploid inds

  # allele frequencies at the loci
  qA <- sum(gtRows[1,])/l
  qB <- sum(gtRows[2,])/l

  qA*qB*(1-qA)*(1-qB)
}

ldDOverSqrtProd <- function(gtRows){
  # takes a 2-row haplotype array
  
  l <- ncol(gtRows) # number of haploid inds
  
  # haplotype frequencies
  xx <- proportions(sapply(0:3, function(x) sum(x == c(1,2) %*% gtRows)))
  
  
  # allele frequencies at the loci
  qA <- sum(gtRows[1,])/l
  qB <- sum(gtRows[2,])/l
  
  
  (xx[4] - qA*qB)/sqrt(qA*qB*(1-qA)*(1-qB))
}

# returns a vector (lower triangle of a matrix)
pairwiseLD <- function(gt){
  #takes haplotype matrix, drops invariant sites
  gts <- gt[apply(gt, 1, sd) != 0,]
  oVals <- matrix(NA, nrow(gts), nrow(gts))
  for(i in 2:nrow(gts)){
    #print(paste0(i, " out of ", nrow(gts)))
    for(j in 1:(i-1)){
      # need to decide which function to call here!
      oVals[i,j] <- ldD(gts[c(i,j),])
    }
  }
  return(oVals[lower.tri(oVals)])
  
}

pairwisePQProd <- function(gt){
  #takes haplotype matrix, drops invariant sites
  gts <- gt[apply(gt, 1, sd) != 0,]
  oVals <- matrix(NA, nrow(gts), nrow(gts))
  for(i in 2:nrow(gts)){
    for(j in 1:(i-1)){
      # need to decide which function to call here!
      oVals[i,j] <- pqProd(gts[c(i,j),])
    }
  }
  oVals[lower.tri(oVals)]
  
}


pairwiser <- function(gt){
  gts <- gt[apply(gt, 1, sd) != 0,]
  as.vector(cor(t(gts))[lower.tri(matrix(0,nrow(gts),nrow(gts)))])
}

meanN <- function(vec, n){
  sum(vec)/n
}

getPBar <- function(gts, l){
  sum(rowSums(gts)/ncol(gts))/l
}

getSitePis <- function(gts){
  p <- rowSums(gts)/ncol(gts)
  return(1 - (p^2 + (1-p)^2))
}

secCol <- function(x){
  a <- strsplit(x, ",")[[1]][[2]]
  n <- nchar(a)
  as.numeric(substr(a, 1, n-1))
}


# neutral expectation (no HRI) as list
sfsExp <- lapply(c(10, 20, 40, 80), function(x){
  1/((1:(x-1))*a_sub_1(x))
})

meanSE <- function(vec){
  m <- mean(vec)
  se <- sd(vec) / sqrt(length(vec))
  c(mean=m, CIlo = m + qnorm(0.025)*se, CIhi = m + qnorm(0.975)*se)
}

a_sub_1 = function(n) {
  sum(1/(1:(n-1)))
}
sum(1/(1:9*a_sub_1(10)))

expNeut <- function(n){
  1/(1:(n-1)*a_sub_1(n))
}
