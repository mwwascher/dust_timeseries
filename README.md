# dust_timeseries

This repository contains the code used to produce figs #### in ####. We use julia to run the MCMC algorithm described in appendix ### and an R script to generate the plots.

The code to run the calibration step is contained in MCMC_timeseries.jl.

The code to run the estimation step is contained in MCMC_timeseries_est.jl.

The R script dust_plots.R contains code to generate the plots from the csv files output from the julia code.
