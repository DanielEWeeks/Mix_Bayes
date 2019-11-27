library(actuar)
library(LearnBayes)
## -----------------------------------
Cover = 0.95
## ----------------------
## rs9469031
## ----------------------
N1 = 1291 + 50
N2 = 1820 + 140 + 2
N = 1 / (1/(2*N1) + 1/(2*N2)) # use allele counts which are twice the numberof subjects
sN = sqrt(N)
p.obt = 1.54e-4
z.obt = qnorm(p.obt/2, lower.tail = F)
snp.maf = (0.019 * N1 + 0.036 * N2) / (N1 + N2)
Sgm = (log(0.73) - log(0.52))/qnorm(1 - 0.05/2) # approx 1/sqrt(N*(snp.maf*(1 - snp.maf)))


## ----------------------------------
## rs200847762
## ----------------------------------
## N1 = 1329 + 12
## N2 = 1150 + 85 + 10
## N = 1 / (1/(2*N1) + 1/(2*N2)) # use allele counts which are twice the numberof subjects
## sN = sqrt(N)
## ## 
## p.obt = 1.84e-6
## ## 
## z.obt = qnorm(p.obt/2, lower.tail = F)
## ##
## snp.maf = (0.024 * N1 + 0.042 * N2) / (N1 + N2)
## ## (Sgm <- sqrt((snp.maf*(1 - snp.maf))^(-1))) ## 5.625543
## Sgm = (log(0.39) - log(0.21))/qnorm(1 - 0.05/2)
## ## Sgm = 1 / Sgm                           #3.166139

dnorm2<-Vectorize(dnorm,vectorize.args="mean")
fu1 <- function(z.obt, fdn, w) {
    w*dnorm2(x=z.obt, mean=fdn) / sum(w*dnorm2(x=z.obt,mean=fdn))
}

## Compute.prior = FALSE
Compute.prior = TRUE
if(Compute.prior) {
    null.pr = 0.99
    null.val = 1.001
    OR = c(1.07, 1.08, 0.92, 0.87, 1.18, 1.181, 1.11, 1.12,
           1.13, 1.12, 1.14, 0.831, 0.86, 1.16, 1.21, 1.20,
           1.19, 0.83, 1.25, 1.26)
    log.OR = log(OR)
    dd = data.frame(log.OR, OR)
    dd <- dd[order(dd$log.OR),]
    nlo = c(0.006006006, 0.006006006, 0.007507508, 0.088588589,
            0.184684685, 0.348348348, 0.205705706, 0.031531532,
            0.016516517, 0.016516517, 0.015015015, 0.015015015,
            0.013513514, 0.012012012, 0.006006006, 0.006006006,
            0.006006006, 0.006006006, 0.004504505, 0.004504505)
    nlo <- nlo/sum(nlo)
    dd$b <- nlo
    dd$ind <- seq_len(nrow(dd))
    dd <- rbind(dd, data.frame(log.OR = log(1/null.val), OR = 1/null.val, b=(null.pr/2), ind=5.1))
    dd <- rbind(dd, data.frame(log.OR = log(null.val), OR = null.val, b=(null.pr/2), ind=5.2))
    dd <- dd[order(dd$ind),]
    v1 = subset(dd, ind<=5 | ind>=6)
    v2 = subset(dd, ind>5 & ind<6)
    v3 = (1 - sum(v2$b)) * v1$b/sum(v1$b)
    v1$b <- v3
    dd[which(dd$ind>=6 | dd$ind<=5),] <- v1
    OR = dd$OR
    log.OR = dd$log.OR
    nlo = dd$nlo
    b = dd$b
    Generate.Prior <- function(w, gammas, oo) {
        Nb = length(gammas)
        idx = sample(1:Nb, prob=w, size=oo, replace=TRUE)
        ## Sg = log(1.02 + 24 / ((1/abs(gammas[idx]))^2))
        ## Sg = 1.75*log(1.002 + 70 / ((1/abs(gammas[idx]))^2))  #<--- Lung_NonConj.RData
        ## Sg <- log( 1 + 9 / ((1/dd$log.OR)^2) ) #<--- Lung_NonConj_1.RData
        ## Sg = 2*log(1 + 70 / ((1/abs(gammas[idx]))^2))  #<--- Lung_NonConj_2.RData
        ## Sg <- log( 1 + 10 / ((1/dd$log.OR)^2) ) #<--- Lung_NonConj_3.RData
        ## Sg <- log( 1 + 9 / ((1/dd$log.OR)^2) ) #<--- Lung_NonConj_4.RData
        Sg <- 0.75 * log( 1 + 9 / ((1/dd$log.OR)^2) ) #<--- Lung_NonConj_5.RData
        zz = rnorm(oo, mean = gammas[idx], sd=Sg)    
        x = zz-mean(zz)
        ncuts = 3000
        cuts = seq(min(x), max(x), length.out = ncuts)
        cat.x = cut(x, cuts)
        lOR.effs = cuts[-1] + (cuts[1:(length(cuts)-1)] - cuts[-1])/2
        lOR.count = as.numeric(table(cat.x))
        lOR.effs = lOR.effs[which(lOR.count>0)]
        lOR.count = lOR.count[which(lOR.count>0)]
        lOR.freq <- lOR.count/sum(lOR.count)
        lOR.effs = lOR.effs - mean(lOR.freq * lOR.effs)
        rm(zz); rm(x)
        rbind(lOR.effs, lOR.freq)
    }
    Prior <- Generate.Prior(b, log.OR, 4e6)
    px <- Prior[2,]
    fud <- Prior[1,]
    save(file="Lung_NonConj_5.RData", px, fud)
} else {
    load(file="Lung_NonConj.RData")
}

b <- px
log.OR <- fud
# sum(px*(fud - mean(fud))^2) # prior variance
sum(px*(fud - sum(px*fud) )^2) # prior variance
par(mfrow=c(1,1)); plot(px ~ fud, xlim=c(-0.2,0.2))
#summary(exp(fud)); summary(1/exp(fud))
x = sample(fud, prob=px, size=3e6, replace=TRUE)
#qu=0.05; ecdf(x)(qu) - ecdf(x)(0); ecdf(x)(0) - ecdf(x)(-qu)
summary(exp(x)); summary(1/exp(x))
quantile(exp(x),0.999); quantile(1/exp(x),0.999)
qu=log(1.01);ecdf(x)(qu)-ecdf(x)(-qu)
qu=log(1.05);ecdf(x)(qu)-ecdf(x)(-qu)
#ecdf(x)(qu) - ecdf(x)(-qu); pnorm(qu, mean(x), sd(x)) - pnorm(-qu, mean(x), sd(x))
# compute noncentralities
z.nonc = log.OR / Sgm # * sqrt(maf*(1-maf) * N)
zz = -z.obt # TAKE NEGATIVE if OR<1
Pst <- fu1(zz, z.nonc, b)
post.mu <- (sum(Pst * z.nonc)) * Sgm #posterior mean for log(OR)
exp(post.mu)
ddist <- cbind(z.nonc, Pst)
hpd.set <- (discint(ddist, Cover)$set) * Sgm
mxbL <- hpd.set[1]                      #posterior LB for log(OR)
mxbU <- hpd.set[length(hpd.set)]        #UB for log(OR)
exp(mxbL)
exp(mxbU)
