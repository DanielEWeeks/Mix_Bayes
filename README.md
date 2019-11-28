* R scripts for sample calculation of
"Flexible Approximate Bayes", FAB

* Chi_square_ln_odds_ratio.r is an illustration of posterior interval
and point estimates for the squared ln(OR). The prior on ln(OR)^2 is
assumed to be the squared normal N(0, s0), that is the scaled
chi-square, s0 \times X, X ~ 1 df chi square or
Gamma(shape = 0.5, scale = 2s0).  L population values of ln(OR)^2
are generated from the prior distribution as 
rgamma(L, shape = 0.5, scale = (2*s0)), binomial samples are taken, 
and L chi-square statistics for testing H0: ln(OR)=0 are computed.
The interval and the point estimates are constructed via FAB for the
largest of the L statistics.

* Wakefield.equivalence.check.r is an illustration of posterior
computation for ln(OR) using FAB and ABF. Normal prior is assumed,
therefore the two methods are equivalent here.

* T-test.r is an illustration of FAB via simulated two-sample multiple
T-tests with selection of the largest T-statistic.

* F-test.r is an illustration of FAB via simulated F-tests, assuming the
prior distribution for mu^2 to follow a Gamma distribution. The
largest of L F-statistics is selected.

* Empirical_prior_plot.r creates an empirical prior distribution for the
analysis of lung cancer data based on frequencies reported in Park et
al. and saves the plot of that distribution into
"Empirical_prior.pdf".
