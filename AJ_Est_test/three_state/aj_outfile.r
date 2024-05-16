args = commandArgs(TRUE)
exact_time = as.logical(as.numeric(args[1]))

par_index = list( beta=1:18)

par_est_mat_split = matrix(nrow = 100, ncol = length(par_index$beta))
par_est_mat       = matrix(nrow = 100, ncol = length(par_index$beta))

# # Load the MCMC results
# chain_list <- NULL
# n_post = 9000; burnin = 1000; steps = 10000
# index_post = (steps - burnin - n_post + 1):(steps - burnin)
# 
# ind_i = 1
# for(i in 1:100) {
#     if(!(i %in% c(39, 61))) {
#         if(exact_time) {
#             load(paste0('Model_out/exactTime/mcmc_out_', i, '.rda'))
#         } else {
#             load(paste0('Model_out/interTime/mcmc_out_', i, '.rda'))
#         }
#         
#         chain_list[[ind_i]] = mcmc_out$chain
#         print(mcmc_out$accept)
#         
#         par_est_mat[ind_i,] = colMeans(mcmc_out$chain[,par_index$beta])   
#         ind_i = ind_i + 1
#     }
# }
# 
# stacked_chains = do.call( rbind, chain_list)

# Load the AJ results
for(i in 1:100) {
    if(exact_time) {
        load(paste0('Model_out/exactTime/par_est_list_', i, '.rda'))    
    } else {
        load(paste0('Model_out/interTime/par_est_list_', i, '.rda'))    
    }
    
    par_est_mat_split[i,] = c(par_est_list[[1]])
}

labels <- c('baseline S1 ---> S2',
            'baseline S1 ---> S3',
            'baseline S2 ---> S1',
            'baseline S2 ---> S3',
            'baseline S3 ---> S1',
            'baseline S3 ---> S2',
            'time S1 ---> S2',
            'time S1 ---> S3',
            'time S2 ---> S1',
            'time S2 ---> S3',
            'time S3 ---> S1',
            'time S3 ---> S2',
            'sex S1 ---> S2',
            'sex S1 ---> S3',
            'sex S2 ---> S1',
            'sex S2 ---> S3',
            'sex S3 ---> S1',
            'sex S3 ---> S2')

trueValues=c(matrix(c(-3,  0.50, -0.4689827,
                      -3,  0.50,  0.2557522,
                      -3,  0.50, -0.1457067,
                      -3,  0.50, -0.8164156,
                      -3,  0.50,  0.5966361,
                      -3,  0.50,  0.7967794), ncol = 3, byrow = T))

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
    yVar = c(par_est_mat_split[,r])#, par_est_mat[,r])
    disc_type = c(rep('AJ', nrow(par_est_mat_split)))#, rep('deSolve', nrow(par_est_mat)))
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

grid.arrange(VP[[1]], VP[[2]], VP[[3]], VP[[4]], VP[[5]], VP[[6]], ncol=2, nrow =3)
grid.arrange(VP[[7]], VP[[8]], VP[[9]], VP[[10]], VP[[11]], VP[[12]], ncol=2, nrow =3)
grid.arrange(VP[[13]], VP[[14]], VP[[15]], VP[[16]], VP[[17]], VP[[18]], ncol=2, nrow =3)

# # Only thinning down the visualization. All means and medians are calculated
# # using the entire stacked_chains matrix.
# thin_ind = seq(1, length(index_post), by=10)
# thin_ind_big = seq(1, nrow(stacked_chains), by=10)
# par_mean = par_median = upper = lower = rep( NA, ncol(stacked_chains))
# par(mfrow=c(4, 2))
# for(r in 1:length(labels)){
# 
#     plot( NULL, xlab=NA, ylab=NA, main=labels[r], xlim=c(1,length(thin_ind)),
#           ylim=range(stacked_chains[,r]) )
# 
#     for(seed in 1:length(chain_list)) lines( chain_list[[seed]][thin_ind,r], type='l', col=seed)
# 
#     par_mean[r] = round( mean(stacked_chains[,r]), 4)
#     par_median[r] = round( median(stacked_chains[,r]), 4)
#     upper[r] = round( quantile( stacked_chains[,r], prob=.975), 4)
#     lower[r] = round( quantile( stacked_chains[,r], prob=.025), 4)
# 
#     print(paste(labels[r], ": [", lower[r], ", ", upper[r], "]"))
# 
#     hist( stacked_chains[thin_ind_big,r], breaks=sqrt(nrow(stacked_chains[thin_ind_big,])), ylab=NA, main=NA,
#           freq=F, xlab=paste0('Mean = ',toString(par_mean[r]),
#                               ' Median = ',toString(par_median[r])))
#     abline( v=upper[r], col='red', lwd=2, lty=2)
#     abline( v=trueValues[r], col='green', lwd=2, lty=2)
#     abline( v=lower[r], col='purple', lwd=2, lty=2)
# }

dev.off()
