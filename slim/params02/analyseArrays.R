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
# gt2 <- read.gt("~/Desktop/sim000.gt2")
# head(gt2)
# dimnames(gt2)
# plot(rownames(gt2))

# are the positions unique?
# nrow(gt2) == length(unique(rownames(gt2)))

# gt2[1:10, 1:10]

gt2sfs <- function(gts){
  # gts is assumed to be 0/1 with the 1st col containing positions
  n <- ncol(gts)
  cts <- rowSums(gts)
  sfs <- sapply(0:n, function(k){
    sum(cts==k)
  })
  names(sfs) <- 0:n
  sfs
}
#gt2sfs(gt2)
#gt2sfs(gt2[,1:20])
#plot(gt2sfs(gt2[,1:20]))


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

t(gt2[1:2,])
table(rowSums(t(gt2[1:2,])))
proportions(table(c(1,2) %*% gt2[1:2,]))

# LD for a 2-row matrix
ldD <- function(gtRows){
  l <- ncol(gtRows) # number of haploid inds
  xx <- proportions(sapply(0:3, function(x) sum(x == c(1,2) %*% gtRows)))
  # allele frequencies at the loci
  pA <- sum(gtRows[1,])/l
  pB <- sum(gtRows[2,])/l
  
  # pAB - pA*pB
  xx[4] - pA*pB
}


#ldD(gt2[1:2,])
#ldD(gt2[1:2,])
#ldD(gt2[7:8,])
#proportions(sapply(0:3, function(x) sum(x==t(gt2[7:8,]) %*% c(1,2))))

pairwiseLD <- function(gts){
  oVals <- matrix(NA, nrow(gts), nrow(gts))
  for(i in 2:nrow(gts)){
    for(j in 1:(i-1)){
      oVals[i,j] <- ldD(gts[c(i,j),])
    }
  }
  oVals[lower.tri(oVals)]
  
}
ld2 <- pairwiseLD(gt2)
#image(ld2)
hist(ld2)


pairwiasr2 <- function(gts){
  as.vector(cor(t(gts))[lower.tri(matrix(0,nrow(gts),nrow(gts)))])
}

#cc <- as.vector(cor(t(gt2))[lower.tri(matrix(0,nrow(gt2),nrow(gt2)))])
#cc2 <- as.vector(cor(t(gt2))[upper.tri(matrix(0,nrow(gt2),nrow(gt2)))])



plot(cc, ld2)
hist(cc)
hist(cc2)
?lower.tri
tajd(gt2sfs(gt2[,1:20]))
theta_pi(gt2sfs(gt2[,1:20]), persite = T)
theta_w(gt2sfs(gt2[,1:20]), persite = T)

deltaTheta(gt2sfs(gt2[,1:20]))


tajd(gt2sfs(gt2[,20:39]))
