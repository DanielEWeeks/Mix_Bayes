remove(list = ls())
load(file="Lung_NonConj_5.RData")
library(actuar)
library(LearnBayes)
library(spatstat)

## -----------------------------------
Cover = 0.95
## ----------------------
## rs9469031
## ----------------------
## N1 = 1291 + 50
## N2 = 1820 + 140 + 2
## N = 1 / (1/(2*N1) + 1/(2*N2)) # use allele counts which are twice the numberof subjects
## sN = sqrt(N)
## p.obt = 1.54e-4
## z.obt = qnorm(p.obt/2, lower.tail = F)
## snp.maf = (0.019 * N1 + 0.036 * N2) / (N1 + N2)
## Sgm = (log(0.73) - log(0.52))/qnorm(1 - 0.05/2) # approx 1/sqrt(N*(snp.maf*(1 - snp.maf)))

## ----------------------
## rs200847762
## ----------------------
## N1 = 1329 + 12
## N2 = 1908 + 74
## N = 1 / (1/(2*N1) + 1/(2*N2)) # use allele counts which are twice the numberof subjects
## sN = sqrt(N)
## p.obt = 1.84e-6
## z.obt = qnorm(p.obt/2, lower.tail = F)
## snp.maf = (0.004 * N1 + 0.019 * N2) / (N1 + N2)
## Sgm = (log(0.39) - log(0.21))/qnorm(1 - 0.05/2) # approx 1/sqrt(N*(snp.maf*(1 - snp.maf)))

## ----------------------
## rs2298090
## ----------------------
## N1 = 1325 + 16
## N2 = 1904 + 78
## N = 1 / (1/(2*N1) + 1/(2*N2)) # use allele counts which are twice the numberof subjects
## sN = sqrt(N)
## p.obt = 6.16e-5
## z.obt = qnorm(p.obt/2, lower.tail = F)
## snp.maf = (0.006 * N1 + 0.020 * N2) / (N1 + N2)
## Sgm = (log(0.56) - log(0.32))/qnorm(1 - 0.05/2) # approx 1/sqrt(N*(snp.maf*(1 - snp.maf)))

## ----------------------
## rs6141383
## ----------------------
N1 = 1277 + 64
N2 = 1934 + 48
N = 1 / (1/(2*N1) + 1/(2*N2)) # use allele counts which are twice the numberof subjects
sN = sqrt(N)
p.obt = 4.80e-4
z.obt = qnorm(p.obt/2, lower.tail = F)
snp.maf = (0.025 * N1 + 0.012 * N2) / (N1 + N2)
Sgm = (log(2.95) - log(2.00))/qnorm(1 - 0.05/2) # approx 1/sqrt(N*(snp.maf*(1 - snp.maf)))

## ---------------------------------------
## Posterior expectation
## ---------------------------------------
dnorm2<-Vectorize(dnorm,vectorize.args="mean")
fu1 <- function(z.obt, fdn, w) {
    w*dnorm2(x=z.obt, mean=fdn) / sum(w*dnorm2(x=z.obt,mean=fdn))
}
b <- px
log.OR <- fud
z.nonc = log.OR / Sgm # * sqrt(maf*(1-maf) * N)
## zz = -z.obt # TAKE NEGATIVE if OR<1
zz = z.obt # TAKE POSITIVE if OR>1
Pst <- fu1(zz, z.nonc, b)
post.mu <- (sum(Pst * z.nonc)) * Sgm #posterior mean for log(OR)
exp(post.mu)
ddist <- cbind(z.nonc, Pst)
hpd.set <- (discint(ddist, Cover)$set) * Sgm
mxbL <- hpd.set[1]                      #posterior LB for log(OR)
mxbU <- hpd.set[length(hpd.set)]        #UB for log(OR)
exp(mxbL)
exp(mxbU)
plog = data.frame(cbind(z.nonc*Sgm, Pst))
names(plog) <- c("Val", "Pr")
## plog$Val.ordered <- plog$Val[order(plog$Val)]
## all.equal(plog$Val, plog$Val.ordered)
## [1] TRUE <--- values are ordered
sum( subset(plog, Val >mxbL & Val < mxbU)$Pr )

## ----------------------------------
## Posterior median + 95% CI
## ----------------------------------
f = ewcdf(plog$Val, weights = plog$Pr)
post.m = quantile(f, c(0.5))
(idx.median <- with(plog, which(Val == post.m)))
i = 1
repeat{
    idx.ub <- idx.median + i
    idx.lb = idx.median 
    if(sum(plog$Pr[idx.lb:idx.ub]) > Cover )  break;
    idx.lb <- idx.median - i
    if(sum(plog$Pr[idx.lb:idx.ub]) > Cover )  break;
    i = i + 1
}
UB = plog$Val[idx.ub]
LB = plog$Val[idx.lb]
sum(subset(plog, Val >= LB & Val<= UB)$Pr )
c(exp(LB), exp(post.m), exp(UB))  

## one step right and left seperatly

i = 1
repeat{
    idx.ub <- idx.median + i
    if(sum(plog$Pr[idx.median:idx.ub]) > Cover / 2) break;
    i = i + 1
    }
UB = plog$Val[idx.ub]                 #CHOOSE TO INCLUDE THE LAST VALUE OR NOT
## 
i = 1
repeat{
    idx.lb <- idx.median - i 
    if(sum(plog$Pr[idx.lb:idx.median]) > Cover / 2) break;
    i = i + 1
}
LB = plog$Val[idx.lb]
sum(subset(plog, Val >= LB & Val<= UB)$Pr )
c(exp(LB), exp(post.m), exp(UB))

## -----------------------------------
## Posterior mean + 95\% CI
## -----------------------------------
(idx.mu <- which.min(abs(plog$Val - post.mu)))
i = 1
repeat{
    idx.ub <- idx.mu 
    idx.lb = idx.mu - i
    if(sum(plog$Pr[idx.lb:idx.ub]) > Cover )  break;
    idx.ub <- idx.mu + i
    if(sum(plog$Pr[idx.lb:idx.ub]) > Cover )  break;
    i = i + 1
}
UB = plog$Val[idx.ub]
LB = plog$Val[idx.lb]
sum(subset(plog, Val >= LB & Val<= UB)$Pr )
c(exp(LB), exp(post.mu), exp(UB))  


## sum(plog$Pr[694:696]) ## 694 696 832
## [1] 0.5604984
## sum(plog$Pr[694:832])
## [1] 0.9218433

## sum(plog$Pr[552:696]) ## 552 696 840
## [1] 0.5660932
## sum(plog$Pr[552:840])
## [1] 0.9513048



