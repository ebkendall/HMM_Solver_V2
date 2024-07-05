library(AalenJohansen)

# Setting the parameter values to the "true" values which generated the data
trueValues= c(-2.31617310,  -1.28756312,  -1.10116400,  -2.52367543,  -2.10384797,
            0.27050001, -11.65470594,  -0.49306415,   0.28862090,   0.22731278,
            -0.39079609,  -0.05894252,  -0.32509646,   0.48631653,   0.99565810,
            -5.28923943,  -0.90870027,  -2.40751854,  -2.44696544,  -6.52252202,
            -6.24090500)
par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)
beta <- matrix(trueValues[par_index$beta], ncol = 3, byrow = F)

par_est_split = vector(mode = 'list', length = 100)
for(it in 1:100) {
    print(it)
    load(paste0("sim_cav_time_inhomog/DataOut/Continuous/cavData", it, ".rda"))
    eid = unique(cavData$ptnum)
    
    # Format the simulated data into the form for AalenJohansen
    cavData1 = cavData[cavData$sex == 1, ]
    cavData0 = cavData[cavData$sex == 0, ]
    
    eid1 = unique(cavData1$ptnum)
    eid0 = unique(cavData0$ptnum)
    
    cavData_aj1 = vector(mode = 'list', length = length(eid1))
    for(i in 1:length(cavData_aj1)) {
        sub_dat = cavData1[cavData1$ptnum == eid1[i], ]
        cavData_aj1[[i]] = list(times = sub_dat$years, states = sub_dat$state)
    }
    
    cavData_aj0 = vector(mode = 'list', length = length(eid0))
    for(i in 1:length(cavData_aj0)) {
        sub_dat = cavData0[cavData0$ptnum == eid0[i], ]
        cavData_aj0[[i]] = list(times = sub_dat$years, states = sub_dat$state)
    }
    
    # Aalen-Johansen estimator (data split according to covariate) -------------
    fit1_split <- aalen_johansen(cavData_aj1)
    fit0_split <- aalen_johansen(cavData_aj0)
    
    v1_split = list()
    v1_split[[1]] <- unlist(lapply(fit1_split$Lambda, FUN = function(L) L[1,2]))
    v1_split[[2]] <- unlist(lapply(fit1_split$Lambda, FUN = function(L) L[1,4]))
    v1_split[[3]] <- unlist(lapply(fit1_split$Lambda, FUN = function(L) L[2,3]))
    v1_split[[4]] <- unlist(lapply(fit1_split$Lambda, FUN = function(L) L[2,4]))
    v1_split[[5]] <- unlist(lapply(fit1_split$Lambda, FUN = function(L) L[3,4]))
    
    v0_split = list()
    v0_split[[1]] <- unlist(lapply(fit0_split$Lambda, FUN = function(L) L[1,2]))
    v0_split[[2]] <- unlist(lapply(fit0_split$Lambda, FUN = function(L) L[1,4]))
    v0_split[[3]] <- unlist(lapply(fit0_split$Lambda, FUN = function(L) L[2,3]))
    v0_split[[4]] <- unlist(lapply(fit0_split$Lambda, FUN = function(L) L[2,4]))
    v0_split[[5]] <- unlist(lapply(fit0_split$Lambda, FUN = function(L) L[3,4]))
    
    v1_t_split <- fit1_split$t
    v0_t_split <- fit0_split$t
    
    
    # Using optim to get parameter estimates -----------------------------------
    min_residuals <- function(data, par) {
        with(data, sum(( (1/par[2]) * (exp(par[1] + par[2]*t + par[3]*x) - exp(par[1] + par[3]*x)) - y)^2))
    }
    
    optim_coeff_split = matrix(0, nrow = 5, ncol = 3)
    for(i in 1:5) {
        
        df0 = data.frame("y" = v0_split[[i]], "t" = v0_t_split, 
                         "x" = rep(0, length(v0_t_split)))
        df1 = data.frame("y" = v1_split[[i]], "t" = v1_t_split, 
                         "x" = rep(1, length(v1_t_split)))
        
        df_tot = rbind(df0, df1)
        
        init_par = beta[i,]
        
        op = optim(par=init_par, fn=min_residuals, data=df_tot, control = list(maxit = 1e5))
        
        if(op$convergence != 0) {print(paste0("split ", i, " no convergence"))}
        
        optim_coeff_split[i, ] = op$par
    }
    par_est_split[[it]] = optim_coeff_split
}

par_est_mat_split = matrix(nrow = 100, ncol = 15)

for(i in 1:100) {
    par_est_mat_split[i,] = c(par_est_split[[i]])
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

# Plot and save the mcmc trace plots and histograms.
library(tidyverse)
library(gridExtra)
library(latex2exp)

pdf('sim_cav_time_inhomog/Plots/aj.pdf')
VP <- vector(mode="list", length = length(labels))
for(r in 1:length(labels)) {
    # Boxplots for par_est_mat
    yVar = c(par_est_mat_split[,r])
    disc_type = rep('AJ', nrow(par_est_mat_split))
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

dev.off()

post_means = par_est_mat_split
save(post_means, file = paste0("sim_cav_time_inhomog/Plots/post_means_aj.rda"))
