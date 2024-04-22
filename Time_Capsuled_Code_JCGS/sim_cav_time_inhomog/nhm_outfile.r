library(tidyverse)
library(gridExtra)

nFrames = 100

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

par_index = list(beta=1:15, misclass=16:19, pi_logit=20:21)

# The true values are set as the posterior means of the thinned last 15,000 steps
# from running the MCMC routine using the numerical ODE solver (seed 10).
# The true values corresponding to the slope coefficient on time are scaled by
# 3 in order to magnify the effect of the different approaches to modelling 
# continuous time HMMs.
load('real_cav_analysis/Model_out/deSolve/mcmc_out_10.rda')
chain = mcmc_out$chain[10000:25001, ]
ind_keep = seq(1, nrow(chain), by=10)
chain = chain[ind_keep, ]
trueValues = colMeans(chain)
trueValues[6:10] = 3 * trueValues[6:10]

year_data = matrix(data=-1, nrow = nFrames, ncol = 21)

cred_set = vector(mode = "list", length = length(trueValues))

for(i in 1:length(cred_set)) {
    cred_set[[i]] = data.frame('lower' = c(-1), 'upper' = c(-1))
}

# -----------------------------------------------------------------------------
# Load Data and Find Credible Sets --------------------------------------------
# -----------------------------------------------------------------------------

row_ind = 1

for (i in 1:nFrames) {
    file_name = paste0("sim_cav_time_inhomog/Model_out/nhm/nhm_out_", i, ".rda")
    if(file.exists(file_name)) {
        load(file_name)
        year_data[row_ind,] = c(gomp_fit$par)

        # Undoing the centering by 5 from nhm_runfile.r
        year_data[row_ind,1:5] = gomp_fit$par[1:5] - 5*gomp_fit$par[6:10]
        
        inv_hess = solve(gomp_fit$hess)
        for(j in 1:length(cred_set)) {
            upper_lim = year_data[row_ind,j] + 1.96 * sqrt(inv_hess[j,j])
            lower_lim = year_data[row_ind,j] - 1.96 * sqrt(inv_hess[j,j])
            cred_set[[j]][row_ind, ] = c(lower_lim, upper_lim)
        }
        
        row_ind = row_ind + 1
    } else {
        print(paste0("File DNE: ", i))
    }
}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Calculating Coverage --------------------------------------------------------
# -----------------------------------------------------------------------------

cov_df = rep(NA,length(trueValues))

for(i in 1:length(trueValues)) {
    val = trueValues[i]
    covrg = mean(cred_set[[i]]$lower <= val & val <= cred_set[[i]]$upper, na.rm=T)
    cov_df[i] = covrg
    print(paste0("Coverage for ", labels[i], " is: ", covrg))
}
# -----------------------------------------------------------------------------

VP <- vector(mode="list", length = length(labels))
for(i in 1:length(trueValues)) {
    y = data.frame(y = year_data[,i], type = rep("nhm", nrow(year_data)))
    
    xlabel = paste0("Coverage: ", round(cov_df[i], digits = 3))
    
    VP[[i]] = ggplot(y, aes(x = type, y = y)) +
        geom_violin(trim=FALSE) +
        geom_boxplot(width=0.1) +
        ggtitle(labels[i]) +
        ylab('') +
        xlab(xlabel) +
        geom_hline(yintercept=trueValues[i], linetype="dashed", color = "red") +
        theme(text = element_text(size = 7))
    
}

pdf("sim_cav_time_inhomog/Plots/nhm.pdf", onefile = T)
grid.arrange(VP[[1]], VP[[2]], VP[[3]], VP[[4]], VP[[5]],
             VP[[6]], VP[[7]], VP[[8]], VP[[9]], ncol=3, nrow =3)
grid.arrange(VP[[10]], VP[[11]], VP[[12]], VP[[13]], VP[[14]],
             VP[[15]], VP[[16]], VP[[17]], VP[[18]], ncol=3, nrow =3)
grid.arrange(VP[[19]], VP[[20]], VP[[21]], ncol=3, nrow =3)
dev.off()

post_means = year_data

save(post_means, file = paste0("sim_cav_time_inhomog/Plots/post_means_nhm.rda"))
save(cov_df, file = paste0("sim_cav_time_inhomog/Plots/cov_df_nhm.rda"))
