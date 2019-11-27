library(actuar)
library(LearnBayes)
library(ggplot2)

null.pr = 0.99
null.val = 1.001
OR = c(1.07, 1.08, 0.92, 0.87, 1.18, 1.181, 1.11, 1.12,
       1.13, 1.12, 1.14, 0.831, 0.86, 1.16, 1.21, 1.20,
       1.19, 0.83, 1.25, 1.26)
( log.OR <- log(OR) )
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
dd <- dd[order(dd$ind),]
oo <- 2e6
Nb = length(dd$log.OR)
idx = sample(1:Nb, prob=dd$b, size=oo, replace=TRUE)
( Sg <- 0.75 * log( 1 + 9 / ((1/dd$log.OR)^2) ) )
zz = rnorm(oo, mean = dd$log.OR[idx], sd=Sg)    
x = zz-mean(zz)
ncuts = 3000
cuts = seq(min(x), max(x), length.out = ncuts)
cat.x = cut(x, cuts)
lOR.effs = cuts[-1] + (cuts[1:(length(cuts)-1)] - cuts[-1])/2 # midpoints, (cuts_j + cuts_{j+1})/2 
lOR.count = as.numeric(table(cat.x))
lOR.effs = lOR.effs[which(lOR.count>0)]
lOR.count = lOR.count[which(lOR.count>0)]
lOR.freq <- lOR.count/sum(lOR.count)
lOR.effs = lOR.effs - mean(lOR.freq * lOR.effs)
ecdf(x)(0.23) - ecdf(x)(-0.19) # ~ 0.77

sum(lOR.freq*(lOR.effs - sum(lOR.freq*lOR.effs) )^2) # prior variance (0.034)
(log(1.27) / qnorm(1 - 0.1))^2 # with probability 0.9, OR lies between 1/1.27 and 1.27


## ---------------------------
## plot
## ---------------------------
tmp <- data.frame(x = lOR.effs,
                  y = lOR.freq)
tmp <- subset(tmp, x > -0.2 & x < 0.25)
    
p <- ggplot()
p + geom_ribbon(aes(x = x, ymax = y, ymin = 0), alpha = 0.5, data = tmp) +
    geom_ribbon(aes(x = lOR.effs, ymax = lOR.freq, ymin = 0), alpha = 0.3, fill = "red") +
    geom_point(aes(x = log.OR, y = b / 140 ), data = dd, colour = "red", size = 3) +
    scale_y_continuous("Density", sec.axis = sec_axis(~.*140, name = "Frequencies reported by Park et al.")) +
    theme_bw() +
    scale_x_continuous("log(OR)")

dev.copy2pdf(file = "Empirical_prior.pdf")
