library(actuar)
library(LearnBayes)

m0 = 0
W = s0 = 0.1
sd.s0 = sqrt(s0)

xN <- round(runif(1, 10, 500))
n1 <- round(0.5 * xN); n2 <- xN - n1
p1 <- runif(1); p2 <- runif(1)
x11 <- rbinom(1, n1, p1); x12 <- n1 - x11
x21 <- rbinom(1, n2, p2); x22 <- n2 - x21
y11 = x11 + 0.5
y12 = x12 + 0.5
y21 = x21 + 0.5
y22 = x22 + 0.5
trunc <- 0.999999; stp = 0.0001
px <- discretize(pnorm(x, mean=m0, sd=sd.s0),
                 method = "upper",
                 from = qnorm(1-trunc, mean=m0, sd=sd.s0),
                 to = qnorm(trunc, mean=m0, sd=sd.s0), step = stp)
fud <- seq(from=qnorm(1-trunc,mean=m0, sd=sd.s0),
           to=(qnorm(trunc,mean=m0, sd=sd.s0)-stp), by=stp)
px <- px/sum(px)
( Nb <- length(px) )
dnorm2<-Vectorize(dnorm,vectorize.args="mean")
fu1 <- function(z.obt, fdn, w) {
    w*dnorm2(x=z.obt, mean=fdn) / sum(w*dnorm2(x=z.obt,mean=fdn))
}
( N = y11+y12+y21+y22 )
c(y11, y12, y21, y22)
sN = sqrt(N)
s00 = s0 * N
V = (1/y11 + 1/y12 + 1/y21 + 1/y22)
R = W / (W+V) # Wakefield's shrinkage factor
theta = log((y11*y22) / (y12*y21))
Sgms = sqrt( N * (1/y11 + 1/y12 + 1/y21 + 1/y22) )
Z = theta / sqrt(V)
delta = Z * Sgms / sN
pmu1 <- delta * s00/(s00 + Sgms^2) # conj post mean in the script
Wak.Post.Mean <- theta * R  # post mean in Wakefield's notation
Pst <- fu1(Z, sN*fud/Sgms, px)
mb.Mean <- sum(Pst * fud) # Mix Bayes posterior mean
cat(pmu1, Wak.Post.Mean, mb.Mean, "\n")
mb.std <- sqrt(sum(Pst * (fud -  mb.Mean)^2)) # Mix Bayes posterior std
sp1 <- sqrt(Sgms^2 * (s00 / (s00 + Sgms^2)) / N ) # conj post std in the script
sp2 <- sqrt(Sgms^2 * s0/(s0 + (Sgms/sN)^2) / N ) # alternative way, conj post var in the script
Wak.Post.Std <- sqrt(V*R) # post std in Wakefield's notation
cat(sp1, sp2, Wak.Post.Std, mb.std, "\n")
