# LPsim

Files to run bootstrapped lasso simulations in a genetic (and polygenetic!) context.
Coming soon: comments, SCAD regression
The first paper from this project is here:(http://dx.doi.org/10.1515/sagmb-2015-0043)
# Simulation Scripts 
- aggregate.R: Rscript to aggregate different results files when called with args by the shell; more important to local implementation than logic of the simulations
- analyze.R: R script to calculate summary stats across aggregated results; again a matter of local implementation
- cor_resolve.R: A basic set of functions to prune correlated variables
- empirical_illustration.R: Functions used in analyzing the real data set described in a paper (to be cited here later, depending on how it fares in peer review)
- percentile_medians.R: helper functions to get medians and MADs of matrices stored in plain text files
- runsims.R: functions to run the simulations--lasso model fitting, bootstrapping, selection of metaparameter nested within bootstrap resamples, output handling in chosan layout--the core of the simulations
