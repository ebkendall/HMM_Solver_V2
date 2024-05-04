# This script file produces trace plots and histograms of the mcmc output files
library(tidyverse)
library(gridExtra)
library(latex2exp)

index_seeds = 1:100

load('DataOut/trueValues.rda')
true_beta = matrix(trueValues[1:15], ncol = 3)

labels <- c('baseline S1 (well)   --->   S2 (mild)',
            'baseline S1 (well)   --->   S4 (dead)',
            'baseline S2 (mild)   --->   S3 (severe)',
            'baseline S2 (mild)   --->   S4 (dead)',
            'baseline S3 (severe)   --->   S4 (dead)',
            'time S1 (well)   --->   S2 (mild)',
            'time S1 (well)   --->   S4 (dead)',
            'time S2 (mild)   --->   S3 (severe)',
            'time S2 (mild)   --->   S4 (dead)',
            'time S3 (severe)   --->   S4 (dead)')
            # 'sex S1 (well)   --->   S2 (mild)',
            # 'sex S1 (well)   --->   S4 (dead)',
            # 'sex S2 (mild)   --->   S3 (severe)',
            # 'sex S2 (mild)   --->   S4 (dead)',
            # 'sex S3 (severe)   --->   S4 (dead)')

par_est_mat = vector(mode = 'list', length = 2)
for(e in 1:2) par_est_mat[[e]] = matrix(ncol = length(labels), nrow = length(index_seeds))


for(i in index_seeds) {
    load(paste0('DataOut/beta_estimates_', i, '.rda'))
    for(e in 1:2) {
        if(length(c(beta_estimates[[e]])) != 10) print(paste0(e, " ", i))
        par_est_mat[[e]][i, ] = c(beta_estimates[[e]])
    }   
}

# Plot and save the mcmc trace plots and histograms.
pdf(paste0('Plots/par_est.pdf'))

for(e in 1:2) {
    VP <- vector(mode="list", length = length(labels))
    for(r in 1:length(labels)) {
        # Adding the boxplots
        yVar = par_est_mat[[e]][,r]
        disc_type = rep(1, nrow(par_est_mat[[e]]))
        y_label = paste0("sex = ", e-1)
        x_label = paste0("Parameter Value: ", round(trueValues[r], 3))
        truth_par = trueValues[r]
        
        if(e == 2) {
            if(r <= 5) {
                x_label = paste0("Parameter Value: ", round(true_beta[r,1] + true_beta[r,3], 3))
                truth_par = true_beta[r,1] + true_beta[r,3]
            }
        }           

        
        plot_df = data.frame(yVar = yVar, disc_type = disc_type)
        VP[[r]] = ggplot(plot_df, aes(x=disc_type, y = yVar)) +
            geom_violin(trim=FALSE) +
            geom_boxplot(width=0.1) +
            ggtitle(labels[r]) +
            ylab(y_label) +
            xlab(x_label) +
            geom_hline(yintercept=truth_par, linetype="dashed", color = "red") +
            theme(text = element_text(size = 7))
    }
    
    grid.arrange(VP[[1]], VP[[2]], VP[[3]], VP[[4]], VP[[5]], ncol=2, nrow =3)
    grid.arrange(VP[[6]], VP[[7]], VP[[8]], VP[[9]], VP[[10]], ncol=2, nrow =3)
}

dev.off()
