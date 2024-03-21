
source("analyseArrays.R")

gt2 <- read.gt(x = "~/Desktop/sim000.gt2")
gt4 <- read.gt(x = "~/Desktop/sim000.gt4")


samp2 <- lapply(c(10, 20, 40, 80), function(n) sampleGt(gt2, n))
samp4 <- lapply(c(10, 20, 40, 80), function(n) sampleGt(gt4, n))

sfs2 <- lapply(samp2, gt2sfs)
sfs4 <- lapply(samp4, gt2sfs)

sapply(sfs2, theta_pi, persite=T)
sapply(sfs2, theta_w, persite=T)
sapply(sfs2, tajd)
sapply(sfs2, deltaTheta)
