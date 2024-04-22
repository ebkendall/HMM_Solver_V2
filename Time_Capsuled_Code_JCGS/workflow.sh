# This workflow file will run everything needed to reproduce the results in 
# the manuscript as well as the supplementary materials. 
# Much of this workflow.sh is done in an embarrassingly parallel fashion. 

####################
# To reproduce the real CAV data results run:

Rscript real_cav_analysis/data_format.r

for seed in {1..10}
do
Rscript real_cav_analysis/mcmc_runfile_deSolve.r $seed
Rscript real_cav_analysis/mcmc_runfile_expm1.r $seed
Rscript real_cav_analysis/mcmc_runfile_expm2.r $seed
Rscript real_cav_analysis/mcmc_runfile_expm3.r $seed
done

# The code above uses 16 threads in parallel.  The for-loop above can be run in an embarrassingly 
# parallel fashion.  The code for each $seed may take about 20 hours to run.

Rscript real_cav_analysis/msm_runfile.r 5
Rscript real_cav_analysis/nhm_runfile.r 5

# The code above calculates the MLE and standard errors instead of doing MCMC.
# The number next to the script corresponds to the seed.

# The following lines produce trace plots and histograms of the MCMC samples.  See the 
# directory real_cav_analysis/Plots for these figures as well as in the supplementary 
# materials. This also produces data corresponding to the 95% credible or 
# confidence sets 
for seed in {1..4}
do
Rscript real_cav_analysis/mcmc_outfile.r $seed
done

Rscript real_cav_analysis/msm_outfile.r 5
Rscript real_cav_analysis/nhm_outfile.r 5
###################


####################
# To reproduce the simulated CAV data results run the following:
# First we generate 100 synthetic datasets
for seed in {1..100}
do
Rscript sim_cav_time_inhomog/cav_simulate.r $seed
done

# Next we run the MCMC routines for each of the five methods
for seed in {1..100}
do
Rscript sim_cav_time_inhomog/mcmc_runfile_deSolve.r $seed
Rscript sim_cav_time_inhomog/nhm_runfile.r $seed
Rscript sim_cav_time_inhomog/mcmc_runfile_expm1.r $seed
Rscript sim_cav_time_inhomog/mcmc_runfile_expm2.r $seed
Rscript sim_cav_time_inhomog/mcmc_runfile_expm3.r $seed
Rscript sim_cav_time_inhomog/msm_runfile.r $seed
done

# The code above uses 16 threads in parallel.  The for-loop above can be run in an embarrassingly 
# parallel fashion.  The code for each $seed may take over 2 days to run.

# The following lines produce trace plots and histograms of the MCMC samples.  See the 
# directory sim_cav_time_inhomog/Plots for these figures as well as in the supplementary 
# materials. This also produces data corresponding to the coverage, 95% credible
# sets, and posterior means.
for seed in {1..4}
do
Rscript sim_cav_time_inhomog/mcmc_outfile.r $seed
done

Rscript sim_cav_time_inhomog/msm_outfile.r
Rscript sim_cav_time_inhomog/nhm_outfile.r
###################


####################
# To reproduce the ECOG sleep analysis results, run the following:

# First we run the MCMC routines for each of the two methods
for seed in {1..100}
do
Rscript real_ecog_analysis/mcmc_runfile.r $seed
Rscript real_ecog_analysis/mcmc_runfile_expm.r $seed
done

# The code above uses 20 threads in parallel.  The for-loop above can be run in an embarrassingly 
# parallel fashion.  The code for each $seed may take over 1.5 days to run.

# The following lines produce trace plots and histograms of the MCMC samples.  See the 
# directory real_ecog_analysis/Plots for these figures. This also produces
# data corresponding to the credible sets.
Rscript real_ecog_analysis/mcmc_outfile.r 1
Rscript real_ecog_analysis/mcmc_outfile.r 2
###################


####################
# To reproduce the computation time results, run the following:
Rscript real_cav_analysis/computation_test.r
Rscript sim_cav_time_inhomog/computation_test.r
Rscript real_ecog_analysis/computation_test.r

# The code above may take 30 minutes to run.
###################

####################
# To reproduce figures 1-5 and all additional figures in the supplementary 
# materials, run the following:
Rscript visualizations/credible_set_visual.r
Rscript visualizations/theorem1_visual.r
Rscript visualizations/violin_plot_visual.r

# To see the plots in the paper, go to the directory visualizations/Plots
# To see the plots in the supplement, go to the directory visualizations/Plots/Supplement
###################
