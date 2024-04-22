# This script file produces trace plots and histograms of the mcmc output files
library(tidyverse)
library(gridExtra)

args = commandArgs(TRUE)
folder = as.numeric(args[1])

model_name = c('deSolve', 'Year', 'YearTwo', 'Month')
dir = paste0('sim_cav_time_inhomog/Model_out/', model_name[folder], '/')

# Size of posterior sample from mcmc chains
n_post = 15000
# Step number at which the adaptive tuning scheme was frozen
burnin = 5000
# Total number of steps the mcmc algorithm is computed for
steps = 20000
# Matrix row indices for the posterior sample to use
index_post = (steps - burnin - n_post + 1):(steps - burnin)

par_index = list(beta=1:15, misclass=16:19, pi_logit=20:21)

index_seeds = 1:100

# The true values are set as the posterior means of the thinned last 15,000 steps
# from running the MCMC routine using the numerical ODE solver (seed 10).
# The true values corresponding to the slope coefficient on time are scaled by
# 3 in order to magnify the effect of the different approaches to modelling 
# continuous time HMMs.
load('real_cav_analysis/Model_out/deSolve/mcmc_out_10.rda')
chain = mcmc_out$chain[10000:25001, ]
ind_keep = seq(1, nrow(chain), by=10)
chain = chain[ind_keep, ]
true_par = colMeans(chain)
true_par[6:10] = 3 * true_par[6:10]

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
# Calculating Credible Sets ---------------------------------------------------
# -----------------------------------------------------------------------------
cred_set = rep(list(matrix(ncol=2,nrow=length(index_seeds))), length(true_par))

ind = 0

for (i in index_seeds) {
    file_name = paste0(dir,'mcmc_out_',toString(i),'.rda')
    if(file.exists(file_name)) {
        load(file_name)
        ind = ind + 1
        for(j in 1:length(true_par)) {

            cred_set[[j]][ind, 1] = round(quantile( mcmc_out$chain[index_post,j],
                                                prob=.025), 4)
            cred_set[[j]][ind, 2] = round(quantile( mcmc_out$chain[index_post,j],
                                                prob=.975), 4)   
        }
    }
}

save(cred_set, file = paste0('sim_cav_time_inhomog/Plots/cred_set_', model_name[folder], '.rda'))

# -----------------------------------------------------------------------------
# Calculating Coverage --------------------------------------------------------
# -----------------------------------------------------------------------------
cov_df = rep(NA, length(true_par))
for(i in 1:length(true_par)) {
    val = true_par[i]
    cov_df[i] = mean(cred_set[[i]][,1] <= val & val <= cred_set[[i]][,2], na.rm=T)
}
save(cov_df, file = paste0('sim_cav_time_inhomog/Plots/cov_df_', model_name[folder], '.rda'))

# -----------------------------------------------------------------------------
# Create mcmc trace plots and histograms
# -----------------------------------------------------------------------------
post_means = matrix(nrow = length(index_seeds), ncol = length(true_par))
chain_list <- NULL
ind = 0
for(seed in index_seeds){
  file_name = paste0(dir,'mcmc_out_',toString(seed), '.rda')
	if(file.exists(file_name)){
	    load(file_name) 
        ind = ind + 1

      	chain_list[[ind]] = mcmc_out$chain[index_post,]
    	post_means[ind,] <- colMeans(mcmc_out$chain[index_post,])
  }
}
save(post_means, file = paste0("sim_cav_time_inhomog/Plots/post_means_",model_name[folder],".rda"))

# Plot and save the mcmc trace plots and histograms.
pdf(paste0('sim_cav_time_inhomog/Plots/mcmc_', model_name[folder], '.pdf'))
par(mfrow=c(4, 2))
stacked_chains = do.call( rbind, chain_list)
par_mean = par_median = upper = lower = rep( NA, ncol(stacked_chains))
VP <- vector(mode="list", length = length(labels))

# Only thinning down the visualization. All means and medians are calculated
# using the entire stacked_chains matrix.
thin_ind = seq(1, length(index_post), by=10)
thin_ind_big = seq(1, nrow(stacked_chains), by=10)

for(r in 1:length(labels)){

	plot( NULL, xlab=NA, ylab=NA, main=labels[r], xlim=c(1,length(thin_ind)),
		    ylim=range(stacked_chains[,r]) )

	for(seed in 1:length(chain_list)) lines( chain_list[[seed]][thin_ind,r], type='l', col=seed)

	par_mean[r] = round( mean(stacked_chains[,r]), 4)
	par_median[r] = round( median(stacked_chains[,r]), 4)
	upper[r] = round( quantile( stacked_chains[,r], prob=.975), 4)
	lower[r] = round( quantile( stacked_chains[,r], prob=.025), 4)
    
    print(paste(labels[r], ": [", lower[r], ", ", upper[r], "]"))

	hist( stacked_chains[thin_ind_big,r], breaks=sqrt(nrow(stacked_chains[thin_ind_big,])), ylab=NA, main=NA,
	      freq=F, xlab=paste0('Mean = ',toString(par_mean[r]),
				                    ' Median = ',toString(par_median[r])))
	abline( v=upper[r], col='red', lwd=2, lty=2)
    abline( v=true_par[r], col='green', lwd=2, lty=2)
	abline( v=lower[r], col='purple', lwd=2, lty=2)
}

for(r in 1:length(labels)) {
    # Adding the boxplots
    yVar = post_means[,r]
    disc_type = rep(model_name[folder], nrow(post_means))
    x_label = paste0("Coverage is: ", round(cov_df[r], digits=3))

    plot_df = data.frame(yVar = yVar, disc_type = disc_type)
    VP[[r]] = ggplot(plot_df, aes(x=disc_type, y = yVar)) +
      geom_violin(trim=FALSE) +
      geom_boxplot(width=0.1) +
      ggtitle(labels[r]) +
      ylab(paste0("Parameter Value: ", round(true_par[r], 3))) +
      xlab(x_label) +
      geom_hline(yintercept=true_par[r], linetype="dashed", color = "red") +
      theme(text = element_text(size = 7))
}

grid.arrange(VP[[1]], VP[[2]], VP[[3]], VP[[4]], VP[[5]],
             VP[[6]], VP[[7]], VP[[8]], VP[[9]], ncol=3, nrow =3)
grid.arrange(VP[[10]], VP[[11]], VP[[12]], VP[[13]], VP[[14]],
             VP[[15]], VP[[16]], VP[[17]], VP[[18]], ncol=3, nrow =3)
grid.arrange(VP[[19]], VP[[20]], VP[[21]], ncol=3, nrow =3)

dev.off()
