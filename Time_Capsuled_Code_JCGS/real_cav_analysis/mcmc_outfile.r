# This script file produces trace plots and histograms of the mcmc output files
library(tidyverse)
library(gridExtra)

args = commandArgs(TRUE)
folder = as.numeric(args[1])

model_name = c('deSolve', 'Year', 'YearTwo', 'Month')
dir = paste0('real_cav_analysis/Model_out/', model_name[folder], '/')

# Size of posterior sample from mcmc chains
n_post = 15000
# Step number at which the adaptive tuning scheme was frozen
burnin = 5000
# Total number of steps the mcmc algorithm is computed for
steps =  30000
# Matrix row indices for the posterior sample to use
index_post = (steps - burnin - n_post + 1):(steps - burnin)

index_seeds = c(1:10) # the number of times the MCMC is run

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)

labels <- c('baseline S1 (well)   --->   S2 (mild)',
            'baseline S1 (well)   --->   S4 (dead)',
            'baseline S2 (mild)   --->   S3 (severe)',
            'baseline S2 (mild)   --->   S4 (dead)',
            'baseline S3 (severe)   --->   S4 (dead)',
            'time S1 (well)   --->   S2 (mild)',
            'time S1 (well)   --->   S4 (dead)',
            'time S2 (mild)   --->   S3 (severe)',
            'time S2 (mild)   --->   S4 (dead)',
            'time S3 (severe)   --->   S4 (dead)',
            'sex S1 (well)   --->   S2 (mild)',
            'sex S1 (well)   --->   S4 (dead)',
            'sex S2 (mild)   --->   S3 (severe)',
            'sex S2 (mild)   --->   S4 (dead)',
            'sex S3 (severe)   --->   S4 (dead)',
            'logit P( obs. state 2 | true state 1 )',
            'logit P( obs. state 1 | true state 2 )',
            'logit P( obs. state 3 | true state 2 )',
            'logit P( obs. state 2 | true state 3 )',
            'logit P( init. state 2 )','logit P( init. state 3 )')

# -----------------------------------------------------------------------------
# Create mcmc trace plots and histograms
# -----------------------------------------------------------------------------
post_means = matrix(nrow = length(index_seeds), ncol = length(labels))
chain_list <- NULL
ind = 0
for(seed in index_seeds){
  file_name = paste0(dir,'mcmc_out_',toString(seed), '.rda')
	if(file.exists(file_name)){
	      load(file_name) 
        ind = ind + 1
        
        # Thinning the chain
        main_chain = mcmc_out$chain[index_post,]
        ind_keep = seq(1, nrow(main_chain), by=10)

      	chain_list[[ind]] = main_chain[ind_keep, ]
    	  post_means[ind,] <- colMeans(main_chain[ind_keep, ])
  }
}

# Plot and save the mcmc trace plots and histograms.
stacked_chains = do.call( rbind, chain_list)
par_mean = par_median = upper = lower = rep( NA, ncol(stacked_chains))
pdf(paste0('real_cav_analysis/Plots/mcmc_', model_name[folder], '.pdf'))
par(mfrow=c(4, 2))

for(r in 1:length(labels)){

	plot( NULL, xlab=NA, ylab=NA, main=labels[r], xlim=c(1,nrow(chain_list[[1]])),
		    ylim=range(stacked_chains[,r]) )

	for(seed in 1:length(chain_list)) lines( chain_list[[seed]][,r], type='l', col=seed)

	par_mean[r] = round( mean(stacked_chains[,r]), 4)
	par_median[r] = round( median(stacked_chains[,r]), 4)
	upper[r] = round( quantile( stacked_chains[,r], prob=.975), 4)
	lower[r] = round( quantile( stacked_chains[,r], prob=.025), 4)
  print(paste(labels[r], ": [", lower[r], ", ", upper[r], "]"))

	hist( stacked_chains[,r], breaks=sqrt(nrow(stacked_chains)), ylab=NA, main=NA,
	      freq=F, xlab=paste0('Mean = ',toString(par_mean[r]),
				                    ' Median = ',toString(par_median[r])))
	abline( v=upper[r], col='red', lwd=2, lty=2)
	abline( v=lower[r], col='purple', lwd=2, lty=2)
}
cred_set_cumulative = cbind(lower, upper)
save(cred_set_cumulative, file = paste0('real_cav_analysis/Plots/cred_set_cumulative_', 
                                            model_name[folder], '.rda'))

dev.off()