# 

library(actuar)
library(LearnBayes)

# location (lo) - scale (s > 0) chi-square transformation
spchisq <- function(x, s, lo=0) { pchisq((x - lo)/s, df=1) } # CDF

sqchisq <- function(x, s, lo=0) { lo + s*qchisq(x, df=1) } # Inverse CDF

sdchisq <- function(x, s, lo=0) { (1/s)*dchisq((x - lo)/s, df=1) } # Density


w.hat = 0.5 # proportion of cases for retrospective sampling
L = 1000 # number of tests
N = 500 # sample size
Cover = 0.95 # target coverage

# parameter for generating L squared ln(OR) from the prior: rgamma(L, shape=0.5, scale=(2*s0))
# (this is equivalent to the scaled 1 df chi-square with the scale factor s0
s0 = 0.3 
sN = sqrt(N)
trunc <- 0.99999

fy.Chi <- function(z, nc) { dchisq(z, df=1, ncp=nc) }
Post.tab.Chi <- function(z, nc, b) {
  nn = length(b)
  dn = fy.Chi(z, nc) %*% b
  pst = rep(0, nn)
  for(i in 1:nn) { pst[i] = fy.Chi(z, nc[i])*b[i] / dn }
  pst
}


## bounds of different frequencies
p.tild.bounds = c(0.1, 0.9) # range of exposure frequency
a1 = 0.9; b1 = 2 # Beta parameters for Pr(D|E), (pDE)
low1 = 0.1; up1 = 0.8 # scale Pr(D|E) to be in this range

Do.Scale <- function(x, x.min, x.max, new.min, new.max) {
    sx = (x - x.min) / (x.max - x.min) # x.min, x.max bounds that x can take
    sx*(new.max - new.min) + new.min # rescale to new bounds
}
## Given a population lnOR, this function generates p1, p2, etc
Generate.Pop.Param <- function(lnOR, w.hat, a1, b1, low1, up1) {
    OR = exp(lnOR)
    k = length(OR)
    pDE <- rbeta(k, a1, b1)
    if(k > 1) {
        pDE <- Do.Scale(pDE, min(pDE), max(pDE), low1, up1) # Pr(D | E)
    }
    pD0 <- 1 / (1 - OR*(1 - 1/pDE)) # Pr(D | non-E)
    p.tild = runif(k, p.tild.bounds[1], p.tild.bounds[2]) # frequency of exposure
    w = pDE*p.tild + pD0*(1-p.tild) # prevalence
    p1 = pDE*p.tild / w
    p2 = (1-pDE)*p.tild / (1-w)
    sigma.pop = sqrt( (1/p.tild) * 1/(pDE*(1-pDE)) + (1/(1-p.tild)) * 1/(pD0*(1-pD0)) )
    sigma.retro = sqrt( (1/w.hat) * 1/(p1*(1-p1)) + (1/(1-w.hat)) * 1/(p2*(1-p2)) )
    delta.pop = lnOR / sigma.pop
    delta.retro = lnOR / sigma.retro
    cbind(p1, p2, w, sigma.pop, sigma.retro, delta.pop, delta.retro, pDE, pD0, p.tild)
}

n1 <- round(w.hat * N); n2 <- N - n1
stp <- 0.0002
px <- discretize(spchisq(x, s=s0), method="round", from=0, to=sqchisq(trunc, s=s0), step=stp)
px <- px/sum(px)

Num.Bins <- length(px)
fud <- seq(from = 0, to = sqchisq(trunc, s=s0) - stp, by = stp)
lOR <- rgamma(L, shape = 0.5, scale = (2*s0))
rv = Generate.Pop.Param(sqrt(lOR), w.hat, a1, b1, low1, up1)
x11 <- rbinom(L, n1, rv[,1]); x12 <- n1 - x11
x21 <- rbinom(L, n2, rv[,2]); x22 <- n2 - x21
y11=x11+0.5; y12=x12+0.5; y21=x21+0.5; y22=x22+0.5
logORs <- ( log((y11*y22) / (y12*y21)) )^2
Sgms <-  ( (N+2)*(1/y11 + 1/y12 + 1/y21 + 1/y22) )
zz <- (N+2) * (logORs / Sgms)
jz <- which.max( zz )
z.obt <- zz[jz]
obs.z <- logORs[jz]
mu.z <- lOR[jz]
fdn2 <- fud * N / Sgms[jz]
pst2 <- Post.tab.Chi(z.obt, fdn2, px)
post.mu.z = (sum(pst2*fdn2)/N) * Sgms[jz]
ddist <- cbind(fud, pst2)
hpd.set <- discint(ddist, Cover)$set
mxbL <- hpd.set[1]
mxbU <- hpd.set[length(hpd.set)]

# Output:
# Lower Interval bound for ln(OR)^2
# Posterior point estimate for ln(OR)^2
# Upper Interval Bound for ln(OR)^2
# True value of ln(OR)^2
cat(mxbL, "<", post.mu.z, "<", mxbU, "; True value:", mu.z, "\n")
