
#setwd("~/git_repos/HRI/slim/params02/")
setwd("~/temp/HRI/")
ff <- dir(pattern = "hri100.o")

an <- function(n){sum(1/1:(n-1))}
an(12)
getNeuSfs <- function(n){
  1/1:(n-1)/an(n)
}
getNeuSfs(12)
readLines(ff[1])

# extract, from SLiM output, the SFS and pBar
sfsAndPbar <- function(f){
  ll <- readLines(f)
  # get last 6 lines (SFSs and stats)
  last6 <- ll[(length(ll)-5):length(ll)]
  
  # split by space and drop 1st element
  l6cl <- lapply(last6, function(x){
    spl <- strsplit(substr(x,2,nchar(x)-1), " ")[[1]]
    as.numeric(spl[2:length(spl)])
  }
  )
  list(SFSsel=l6cl[[1]],
       SFSneu=l6cl[[2]],
       pBarSel=l6cl[[3]],
       pBarNeu=l6cl[[4]],
       PiSel=l6cl[[5]],
       PiNeu=l6cl[[6]]
       )
}

from100 <- as.data.frame(t(sapply(ff, sfsAndPbar)))
head(from100)
from100$pBarSel <- unlist(from100$pBarSel)
from100$pBarNeu <- unlist(from100$pBarNeu)
from100$PiSel <- unlist(from100$PiSel)
from100$PiNeu <- unlist(from100$PiNeu)
hist(from100$pBarSel, breaks="FD")
hist(from100$PiNeu, breaks="FD")
hist(from100$PiSel, breaks="FD", add =T, col="#FF000040")

t.test(from100$PiNeu, from100$PiSel)

coef(summary(lm(from100$pBarSel~1)))
coef(summary(lm(from100$PiNeu~1)))
piNeuHat <- coef(summary(lm(from100$PiNeu~1)))[1,1]
coef(summary(lm(from100$PiSel~1)))
sfsSel <- do.call(rbind, from100$SFSsel)
sfsNeu <- do.call(rbind, from100$SFSneu)
barplot(colMeans(sfsSel)[2:12])
barplot(colMeans(sfsNeu)[2:12])

sum(c(0.413224, 0.177479, 0.105818, 0.0728697, 0.0545042, 0.0430248, 
  0.035276, 0.0297467, 0.0256311, 0.0224646, 0.0199623) * 1:11 * 11:1)/12/12/2500



barplot(rbind(colMeans(sfsSel)[2:12]/sum(colMeans(sfsSel)[2:12]),
      colMeans(sfsNeu)[2:12]/sum(colMeans(sfsNeu)[2:12]),
      piNeuHat/1:11/sum(piNeuHat/1:11),
      c(0.413224, 0.177479, 0.105818, 0.0728697, 0.0545042, 0.0430248, 
        0.035276, 0.0297467, 0.0256311, 0.0224646, 0.0199623)),
beside = T,
names.arg = 1:11,
col=c("#FF000030","#00000040", "#00000080", "#000000"),
main="Avg. SFS (100 replicate simulations)")

abline(h=0:8/20, col = "lightgrey", lty=3)
#grid(nx = NA, ny=8)
legend("top",
       pch=15,
       col=c("#FF000030","#00000040", "#00000080", "#000000"),
       legend=c("Deleterious observed", "Neutral observed", "Neutral (obs. theta/i)", "Neutral expected (HRI notes/P&K)")
       )


# Tajima's D --------------------------------------------------------------
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

tajd(from100$SFSneu[[2]])

from100$DTSel <- sapply(from100$SFSsel, tajd)
from100$DTNeu <- sapply(from100$SFSneu, tajd)

from100$piS <- sapply(from100$SFSsel, theta_pi)/2500
from100$piN <- sapply(from100$SFSneu, theta_pi)/2500

from100$twS <- sapply(from100$SFSsel, theta_w)/2500
from100$twN <- sapply(from100$SFSneu, theta_w)/2500

from100$deltheS <- sapply(from100$SFSsel, function(x) 1- theta_pi(x)/theta_w(x))
from100$deltheN <- sapply(from100$SFSneu, function(x) 1- theta_pi(x)/theta_w(x))



summary(lm(from100$DTSel~1))
summary(lm(from100$DTNeu~1))

summary(lm(from100$piS~1))
summary(lm(from100$piN~1))

summary(lm(from100$twS~1))
summary(lm(from100$twN~1))

summary(lm(from100$deltheN~1))

cbind(sapply(from100$SFSsel, function(x) theta_pi(x)/2500, USE.NAMES = F),
      from100$PiSel)
