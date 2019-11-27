# Illustration of FAB via simulated F-tests, assuming the prior distribution for mu^2 
# to be a Gamma. The largest of L F-statistics is selected

library(actuar)
library(LearnBayes)

L <- 5000 # number of F-tests
N1 <- 20 # sample size, group 1
N2 <- 25 # sample size, group 2
Cover = 0.95 # target coverage
Scale <- 0.5 # prior variance for the distribution of the true difference^2 (mu^2)
Shape <- 1 
trunc <- 0.999999
# Two-sample F test
Sim.Ftest <- function(mu, Sigma, n1, n2) {
  y1 = rnorm(n1, mean = sqrt(mu), sd = Sigma)
  y2 = rnorm(n2, mean = 0, sd = Sigma)
  vr = (((n1-1)*var(y1) + (n2-1)*var(y2))/(n1+n2-2)) # assume equal var
  svr = sqrt(vr)
  den = svr * sqrt(1/n1 + 1/n2)
  mean.diff = (mean(y1)-mean(y2))
  t.stat =  mean.diff / den
  c(t.stat^2, svr^2, mean.diff^2)
}
N = 1/(1/N1 + 1/N2); sN = sqrt(N) # sample size used to compute each F-statistic
df.t = N1+N2-2
#
fy.F <- function(z, df.t, nc) { df(z, df1=1, df2=df.t, ncp = nc) }
Post.tab.F <- function(z, df.t, nc, b) {
  nn = length(b)
  dn = fy.F(z, df.t, nc) %*% b
  pst = rep(0, nn)
  for(i in 1:nn) { pst[i] = fy.F(z, df.t, nc[i])*b[i] / dn }
  pst
}
stp <- 0.001
px <- discretize(pgamma(x, shape=Shape, scale=Scale),
                 method = "round",
                 from = 0,
                 to = qgamma(trunc, shape=Shape, scale=Scale), step = stp)
px <- px/sum(px)
Num.Bins <- length(px)
cat(Num.Bins, "\n")
Num.Bins <- length(px)
fud <- seq(from=0, to=qgamma(trunc, shape=Shape, scale=Scale)-stp, by=stp)
Mu <- rgamma(L, shape=Shape, scale=Scale)
Sigma <- runif(L, Scale/2, Scale+Scale/2)
MoS <- Mu / Sigma
tmp.stat.mx <- c(-Inf, 1, 0)
jmx <- 1
for(j in 1:L) {
    tmp.stat <- Sim.Ftest(Mu[j], Sigma[j], N1, N2)
    if(tmp.stat[1] > tmp.stat.mx[1]) { tmp.stat.mx <- tmp.stat; jmx <- j }
}
t.stat <- tmp.stat.mx[1]
sigma.stat <- tmp.stat.mx[2]
mean.diff <- tmp.stat.mx[3]
Mu.over.Sigma <- MoS[jmx]
Sigma.mx <-  Sigma[jmx]
selected.mu <- Mu[jmx]
fdn <-  (fud*N) / sigma.stat
pst = Post.tab.F(t.stat, df.t, fdn, px)
th.e <- (sum(pst*fdn)/N) * sigma.stat # posterior mean, sigma estimated
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
