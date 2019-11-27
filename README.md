Wakefield.equivalence.check.r is an illustration of posterior
computation for ln(OR) using FAB and ABF. Normal prior is assumed,
therefore the two methods are equivalent here.

T-test.r is an ilustration of FAB via simulated two-sample multiple
T-tests with selection of the largest T-statistic.

F-test.r is an illustration of FAB via simulated F-tests, assuming the
prior distribution for mu^2 to follow a Gamma distribution. The
largest of L F-statistics is selected.

Empirical_prior_plot.r creates an empirical prior distribution for the
analysis of lung cancer data based on frequencies reported in Park et
al. and saves the plot of that distribution into
"Empirical_prior.pdf".
