args = commandArgs(TRUE)
exact_time = as.logical(as.numeric(args[1]))

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)

par_est_mat_split = matrix(nrow = 100, ncol = 15)
par_est_mat       = matrix(nrow = 100, ncol = 15)

# Load the MCMC results
chain_list <- NULL
index_post = 1:4501

ind = 0
for(i in 1:nrow(par_est_mat)) {
    if(exact_time) {
        file_name = paste0('Model_out/exactTime/mcmc_out_', i, '.rda')
    } else {
        file_name = paste0('Model_out/interTime/mcmc_out_', i, '.rda')
    }

    if(file.exists(file_name)) {
        load(file_name)
        mcmc_out$chain[, 1:5] = mcmc_out$chain[, 1:5] - mcmc_out$chain[, 6:10] * mcmc_out$mean_t

        chain_list[[ind+1]] = mcmc_out$chain
    
        par_est_mat[ind+1,] = colMeans(mcmc_out$chain[,par_index$beta])
        ind = ind + 1
    }
}

par_est_mat = par_est_mat[1:ind, ]

stacked_chains = do.call( rbind, chain_list)

# Load the AJ results
for(i in 1:nrow(par_est_mat_split)) {
    if(exact_time) {
        load(paste0('Model_out/exactTime/par_est_list_', i, '.rda'))
    } else {
        load(paste0('Model_out/interTime/par_est_list_', i, '.rda'))
    }

    par_est_mat_split[i,] = c(par_est_list[[1]])
}

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
            'sex S3 (severe)   --->   S4 (dead)')

load('mcmc_out_10.rda')
chain = mcmc_out$chain[10000:25001, ]
ind_keep = seq(1, nrow(chain), by=10)
chain = chain[ind_keep, ]
trueValues = colMeans(chain)
trueValues[6:10] = 3 * trueValues[6:10]

# Plot and save the mcmc trace plots and histograms.
library(tidyverse)
library(gridExtra)
library(latex2exp)

pdf_title = NULL
if(exact_time) {
    pdf_title = 'Plots/par_est_optim_exact.pdf'
} else {
    pdf_title = 'Plots/par_est_optim_inter.pdf'
}

pdf(pdf_title)
VP <- vector(mode="list", length = length(labels))
for(r in 1:length(labels)) {
    # Boxplots for par_est_mat
    yVar = c(par_est_mat_split[,r], par_est_mat[,r])
    disc_type = c(rep('AJ', nrow(par_est_mat_split)), rep('MCMC', nrow(par_est_mat)))
    x_label = paste0("Parameter Value: ", round(trueValues[r], 3))
    truth_par = trueValues[r]
    
    plot_df = data.frame(yVar = yVar, disc_type = disc_type)
    VP[[r]] = ggplot(plot_df, aes(x=disc_type, y = yVar)) +
        geom_violin(trim=FALSE) +
        geom_boxplot(width=0.1) +
        ggtitle(labels[r]) +
        ylab(" ") +
        xlab(x_label) +
        geom_hline(yintercept=truth_par, linetype="dashed", color = "red") +
        theme(text = element_text(size = 7))
}

grid.arrange(VP[[1]], VP[[2]], VP[[3]], VP[[4]], VP[[5]], ncol=2, nrow =3)
grid.arrange(VP[[6]], VP[[7]], VP[[8]], VP[[9]], VP[[10]], ncol=2, nrow =3)
grid.arrange(VP[[11]], VP[[12]], VP[[13]], VP[[14]], VP[[15]], ncol=2, nrow =3)

# Only thinning down the visualization. All means and medians are calculated
# using the entire stacked_chains matrix.
thin_ind = seq(1, length(index_post), by=10)
thin_ind_big = seq(1, nrow(stacked_chains), by=10)
par_mean = par_median = upper = lower = rep( NA, ncol(stacked_chains))
par(mfrow=c(4, 2))
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
    abline( v=trueValues[r], col='green', lwd=2, lty=2)
    abline( v=lower[r], col='purple', lwd=2, lty=2)
}

dev.off()