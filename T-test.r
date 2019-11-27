# Ilustration of FAB via simulated two-sample multiple T-tests with selection of the largest T-statisticv3: changed Sigma <

library(actuar)
library(LearnBayes)

L <- 5000 # number of T-tests
N1 <- 25 # sample size, group 1
N2 <- 20 # sample size, group 2
Cover = 0.95 # target coverage
s0 <- 1 # prior variance for the distribution of the true difference (mu)
trunc <- 0.99999
Sim.Ttest <- function(mu, Sigma, n1, n2) { # two-sample T-test
  y1 = rnorm(n1, mean = mu, sd = Sigma)
  y2 = rnorm(n2, mean = 0, sd = Sigma)
  vr = (((n1-1)*var(y1) + (n2-1)*var(y2))/(n1+n2-2)) # assume equal var
  svr = sqrt(vr)
  den = svr * sqrt(1/n1 + 1/n2)
  mean.diff = (mean(y1)-mean(y2))
  t.stat =  mean.diff / den
  c(t.stat, svr, mean.diff)
}
#
N = 1/(1/N1 + 1/N2); sN = sqrt(N) # sample size used to compute each T
df.t = N1+N2-2
m0 = 0 # prior for mu ~ N(m0, s0)
sd.s0 = sqrt(s0)
#
fy.T <- function(z, df.t, nc) { dt(-z, df=df.t, ncp = (-nc)) }
Post.tab.T <- function(z, df.t, nc, b) {
  nn = length(b)
  dn = fy.T(z, df.t, nc) %*% b
  pst = rep(0, nn)
  for(i in 1:nn) { pst[i] = fy.T(z, df.t, nc[i])*b[i] / dn }
  pst
}
stp <- 0.002
px <- discretize(pnorm(x, mean=m0, sd=sd.s0),
                 method = "round",
                 from = qnorm(1-trunc, mean=m0, sd=sd.s0),
                 to = qnorm(trunc, mean=m0, sd=sd.s0), step = stp)
px <- px/sum(px)
Num.Bins <- length(px)
cat(Num.Bins,"\n")
fud <- seq(from=qnorm(1-trunc,mean=m0, sd=sd.s0), to=qnorm(trunc, mean=m0, sd=sd.s0)-stp, by=stp)
Mu <- rnorm(L, mean=m0, sd=sd.s0)
Sigma <- runif(L, s0-s0/2, s0+s0/2)
MoS <- Mu / Sigma
tmp.stat.mx <- c(-Inf, 1, 0)
jmx <- 1
for(j in 1:L) {
    tmp.stat <- Sim.Ttest(Mu[j], Sigma[j], N1, N2)
    if(tmp.stat[1] > tmp.stat.mx[1]) { tmp.stat.mx <- tmp.stat; jmx <- j }
}
t.stat <- tmp.stat.mx[1]
sigma.stat <- tmp.stat.mx[2]
mean.diff <- tmp.stat.mx[3]
Mu.over.Sigma <- MoS[jmx]
Sigma.mx <-  Sigma[jmx]
selected.mu <- Mu[jmx]
fdn <-  (fud*sN) / sigma.stat
pst <- Post.tab.T(t.stat, df.t, fdn, px) #Post.t(t.stat, df.t, fdn, px)
th.e <- (sum(pst*fdn)/sN) * sigma.stat # posterior mean, sigma estimated
ddist <- cbind(fud, pst)
hpd.set <- discint(ddist, Cover)$set
mxbL <- hpd.set[1]
mxbU <- hpd.set[length(hpd.set)]
# Output:
# Lower Interval bound
# Posterior point estimate for mu
# Upper Interval Bound
# True mean difference
cat(mxbL, "<", th.e, "<", mxbU, "; True value:", selected.mu, "\n")
