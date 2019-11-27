# Illustration of FAB via simulated F-tests, assuming the prior distribution for mu^2 
# to be a Gamma. The largest of L F-statistics is selected to illustrate FAB's resistance
# to selection bias. Possible command line:
# R CMD BATCH --vanilla --slave --no-timing '--args  s0=0.1 B=5000 L=500 N1=100 N2=100' F-test_num_bins_estimates.r /dev/null &
library(actuar)
library(LearnBayes)

# default values
L <- 500 # number of F-tests
N1 <- 25 # sample size, group 1
N2 <- 25 # sample size, group 2
B <- 5000
Cover = 0.80; alpha = 1-Cover
s0 <- 1 # prior variance for the distribution of the true difference^2 (mu^2)
Shape <- 0.5
n.steps <- 25
trunc <- 0.99999
args=(commandArgs(TRUE))
if(length(args) > 0) {
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
#
Sim.Ftest <- function(mu, Sigma, n1, n2, tails=1) {
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
m0 = 0 # prior for mu ~ N(m0, s0)
sd.s0 = sqrt(s0)
#
mean.est <- zz <- th.e <- selected.mu <- rep(0, B)

fy.F <- function(z, df.t, nc) { df(z, df1=1, df2=df.t, ncp = nc) }
Post.tab.F <- function(z, df.t, nc, b) {
  nn = length(b)
  dn = fy.F(z, df.t, nc) %*% b
  pst = rep(0, nn)
  for(i in 1:nn) { pst[i] = fy.F(z, df.t, nc[i])*b[i] / dn }
  pst
}
dnonc2<-Vectorize(df, vectorize.args="ncp")
# Post.F is the same as Post.tab.F but possibly faster
Post.F <- function(z.obt, dfr, fdn, w) {
    w*dnonc2(x = z.obt, df1=1, df2=dfr, ncp = fdn) / sum(w*dnonc2(x=z.obt, df1=1, df2=dfr, ncp=fdn))
}

steps <- rep(0, n.steps)
steps[1] <- 0.0018
for(i in 2 : length(steps)) {
    steps[i] <- steps[i-1] * 1.34
}
steps
#steps <- c(0.000784, 0.0011756, 0.00235, 0.0047, 0.0094, 0.115)

for(jj in steps) {
    stp <- jj
    px <- discretize(pgamma(x, shape=Shape, scale=s0),
                 method = "round",
                 from = 0,
                 to = qgamma(trunc, shape=Shape, scale=s0), step = stp)
    px <- px/sum(px)
    Num.Bins <- length(px)
    cat(Num.Bins, "")
}
cat("\n")

File1 <- paste("F_test_est_Num.Bins.", "s0.", s0, ".N1.", N1, ".N2.", N2, ".L.", L , ".txt", sep="")
cat(file=File1, "Num.Bins,", "True mean,", "Posterior mean", "\n", sep="")

for(jj in steps) {
    stp <- jj
    px <- discretize(pgamma(x, shape=Shape, scale=s0),
                     method = "round",
                     from = 0,to = qgamma(trunc, shape=Shape, scale=s0), step = stp)
    px <- px/sum(px)
    Num.Bins <- length(px)
    fud <- seq(from=0, to=qgamma(trunc, shape=Shape, scale=s0)-stp, by=stp)
    mxb <- 0
    i <- 1
    repeat {
    	if(i > B) break;
        Mu <- rgamma(L, shape=Shape, scale=s0)
        Sigma <- runif(L, s0/2, s0+s0/2)
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
        zz[i] <- t.stat
        selected.mu[i] <- Mu[jmx]
        fdn <-  (fud*N) / sigma.stat # <-- IMPORTANT: when prior is on Mu, need to recompute "fdn" (noncentrality) each time
        pst = Post.tab.F(t.stat, df.t, fdn, px)
        th.e[i] <- (sum(pst*fdn)/N) * sigma.stat # posterior mean, sigma estimated
        mean.est[i] <- mean.diff
        if(!is.na(th.e[i]) && !is.na(selected.mu[i])) {
          #ddist <- cbind(fud, pst)
          #hpd.set <- discint(ddist, Cover)$set
          #mxbL <- hpd.set[1]
          #mxbU <- hpd.set[length(hpd.set)]
          #if(Mu[jmx] > mxbL && Mu[jmx] < mxbU) mxb <- mxb+1
          if(i %% 10 == 0) cat(".", i, ".", sep=""); if(i %% 50 == 0) cat("\r")
          i <- i+1
        }
    }
    cat("\n")
    #emp.Cov <- sum(mxb)/B
    #cat(file=File1, Num.Bins, ",", stp, ",",  emp.Cov, "\n", sep="", append=TRUE)
    cat(file=File1, Num.Bins, ",", stp, ",",  mean(selected.mu), ",", mean(th.e), "\n", sep="", append=TRUE)
    #cat(Num.Bins, ",", stp, ",",  mean(selected.mu), ",", mean(th.e), "\n", sep="")
}

#read.csv(File1, as.is=TRUE)
