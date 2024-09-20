

#gg <- read.gt("params02/data/sim0197.gt2")
gg <- read.gt("params03/data/sim0200.gt2")
gg[1:10, 1:10]
ss <- 1e-2 # selecton coef (although for deleterious variants, it is positive in the model!)
nn <- 2500 # number of sites in the genome

wSum <- colSums(gg * -ss ) + 1
wProd <- exp(colSums(log(gg * -ss +1)))
var(wSum)
var(wProd)
var(log(wProd))
hist(wSum)
hist(wProd)
var(wSum)
var(wProd)
# fitness variance per site
siteVars <- apply(gg * -ss + 1, 1, var)
# the sum of those
vg <- sum(siteVars)

# total var in fitnes across individuals
va <- var(colSums(gg * -ss))
va / vg
sum(getSitePis(gg)) * ss * ss /2 # expectd var in fitness
vg
va
(va - vg) / (2 * ss * ss)
ld <- pairwiseLD(gg)
sum(ld) / nn *2

ldD(gg[6:7,])
gg[1:10, 1:10]
var(gg[1,])
var(gg[2,])
gg[6:7,]

