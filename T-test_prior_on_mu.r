# Ilustration of FAB via simulated two-sample multiple T-tests with selection of the largest T-statistic

# Possible command line supplying script parameters:
# R CMD BATCH --vanilla --slave --no-timing '--args  s0=0.5 N1=20 N2=30 mx=1000' T-test_prior_on_mu.r /dev/null &
library(actuar)
library(LearnBayes)

# default values
mx <- 2500 # number of T-tests
N1 <- 50 # sample size, group 1
N2 <- 70 # sample size, group 2
B <- 1000
s0 <- 0.1 # prior variance for the distribution of the true difference (mu)
args=(commandArgs(TRUE))
if(length(args) > 0) {
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
# two-sample T-test
Sim.Ttest <- function(mu, Sigma, n1, n2, tails=1) {
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
N = 1/(1/N1 + 1/N2); sN = sqrt(N) # sample size used to compute each T-statistic
df.t = N1+N2-2
m0 = 0 # prior for mu ~ N(m0, s0)
sd.s0 = sqrt(s0)
#
mean.est <- zz <- th.t <-th.e <- selected.mu <- rep(0, B)
trunc <- 0.9999; stp = 0.01
# assume normal prior for mu
px <- discretize(pnorm(x, mean=m0, sd=sd.s0),
                 method = "upper",
                 from = qnorm(1-trunc, mean=m0, sd=sd.s0),
                 to = qnorm(trunc, mean=m0, sd=sd.s0), step = stp)
px <- px/sum(px)
fud <- seq(from=qnorm(1-trunc,mean=m0, sd=sd.s0), to=qnorm(trunc, mean=m0, sd=sd.s0)-stp, by=stp)
#fdn = sN * fud # <-- THIS IS OK WHEN THE PRIOR IS ON MU/SIGMA
#
fy.T <- function(z, df.t, nc) { dt(z, df=df.t, ncp=nc) }
Post.tab.T <- function(z, df.t, nc, b) {
  nn = length(b)
  dn = fy.T(z, df.t, nc) %*% b
  pst = rep(0, nn)
  for(i in 1:nn) { pst[i] = fy.T(z, df.t, nc[i])*b[i] / dn }
  pst
}
dnonc2<-Vectorize(dt,vectorize.args="ncp")
# Post.t is the same as Post.tab.T but possibly faster
Post.t <- function(z.obt, dfr, fdn, w) {
    w*dnonc2(x=z.obt, df=dfr, ncp=fdn) / sum(w*dnonc2(x=z.obt, df=dfr, ncp=fdn))
}
#
i <- 1
for(i in 1:B) {
    Mu <- rnorm(mx, mean=m0, sd=sd.s0)
    Sigma <- runif(mx, 1, 5)
    MoS <- Mu / Sigma
    tmp.stat.mx <- c(-Inf, 1, 0)
    jmx <- 1
    for(j in 1:mx) {
        tmp.stat <- Sim.Ttest(Mu[j], Sigma[j], N1, N2)
        if(tmp.stat[1] > tmp.stat.mx[1]) { tmp.stat.mx <- tmp.stat; jmx <- j }
    }
    t.stat <- tmp.stat.mx[1]
    sigma.stat <- tmp.stat.mx[2]
    mean.diff <- tmp.stat.mx[3]
    Mu.over.Sigma <- MoS[jmx]
    Sigma.mx <-  Sigma[jmx]
    zz[i] <- t.stat
    selected.mu[i] <- Mu[jmx]
    fdn <-  (fud*sN) / sigma.stat # <-- IMPORTANT: when prior is on Mu, need to recompute "fdn" (noncentrality) each time
    pst = Post.t(t.stat, df.t, fdn, px)
    th.t[i] <- (sum(pst*fdn)/sN) * Sigma.mx # posterior mean, sigma known
    th.e[i] <- (sum(pst*fdn)/sN) * sigma.stat # posterior mean, sigma estimated
    mean.est[i] <- mean.diff
    if(i %% 10 == 0) cat(".", i, ".", sep=""); if(i %% 50 == 0) cat("\r")
}
cat("\n")
#OutFile = paste("T-test.max.mu.", "mx_", mx, ".s0_", s0, ".N_", N, ".txt", sep="")
OutFile = "/dev/stdout"
cat(file=OutFile, "num bins =", length(px), "\n", append=TRUE)
cat(file=OutFile, "Number of tests (L) =", mx, "\n", append=TRUE)
#cat(file=OutFile, "Bayes, sigma known:", 100*Count.D/B,  "Bayes, sigma estimated:", 100*Count.E/B, append=TRUE)
cat(file=OutFile, "average post, s. known =", mean(th.t,na.rm=TRUE), "\n",
    "average post, s. estimated: =", mean(th.e,na.rm=TRUE), "\n",
    "average frequentist estimated: =", mean(mean.est,na.rm=TRUE), "\n",
    "actual =", mean(selected.mu), "\n", append=TRUE)
par(mfrow=c(2,2))
plot(th.e, selected.mu)
plot(th.t, selected.mu)
plot(mean.est, selected.mu)
plot(th.t, th.e)

#shft = 0.5 + abs(min(c(th.e, th.t, selected.mu, mean.est)))
#plot(log(shft+th.e), log(shft+selected.mu))
#plot(log(shft+th.t), log(shft+selected.mu))
#plot(log(shft+mean.est), log(shft+selected.mu))
#plot(log(shft+th.t), log(shft+th.e))
