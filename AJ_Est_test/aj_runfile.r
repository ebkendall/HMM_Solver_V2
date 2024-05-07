library(AalenJohansen)

args = commandArgs(TRUE)
exact_time = as.logical(as.numeric(args[1]))

load(paste0('DataOut/trueValues.rda'))
par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)
beta <- matrix(trueValues[par_index$beta], ncol = 3, byrow = F)

# Format the simulated data into the form for AalenJohansen
par_est = vector(mode = 'list', length = 100)
par_est_split = vector(mode = 'list', length = 100)
for(it in 1:100) {
    print(it)
    if(exact_time) {
        load(paste0('DataOut/exactTime/cavData', it, '.rda'))
    } else {
        load(paste0('DataOut/interTime/cavData', it, '.rda'))
    }
    eid = unique(cavData$ptnum)
    cavData_aj = vector(mode = 'list', length = length(eid))
    for(i in 1:length(cavData_aj)) {
        sub_dat = cavData[cavData$ptnum == eid[i], ]
        cavData_aj[[i]] = list(times = sub_dat$years, states = sub_dat$state, X = sub_dat$sex[1])
    }
    
    
    # Conditional Aalen-Johansen estimator ------------------------------------
    fit1 <- aalen_johansen(cavData_aj, x = 1)
    fit0 <- aalen_johansen(cavData_aj, x = 0)
    
    v1 = list()
    v1[[1]] <- unlist(lapply(fit1$Lambda, FUN = function(L) L[1,2]))
    v1[[2]] <- unlist(lapply(fit1$Lambda, FUN = function(L) L[1,4]))
    v1[[3]] <- unlist(lapply(fit1$Lambda, FUN = function(L) L[2,3]))
    v1[[4]] <- unlist(lapply(fit1$Lambda, FUN = function(L) L[2,4]))
    v1[[5]] <- unlist(lapply(fit1$Lambda, FUN = function(L) L[3,4]))
    
    v0 = list()
    v0[[1]] <- unlist(lapply(fit0$Lambda, FUN = function(L) L[1,2]))
    v0[[2]] <- unlist(lapply(fit0$Lambda, FUN = function(L) L[1,4]))
    v0[[3]] <- unlist(lapply(fit0$Lambda, FUN = function(L) L[2,3]))
    v0[[4]] <- unlist(lapply(fit0$Lambda, FUN = function(L) L[2,4]))
    v0[[5]] <- unlist(lapply(fit0$Lambda, FUN = function(L) L[3,4]))
    
    v1_t <- fit1$t
    v0_t <- fit0$t
    
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
    
    # Unconditional Aalen-Johansen estimator (but data is split) ---------------
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
    
    optim_coeff = matrix(0, nrow = 5, ncol = 3)
    for(i in 1:5) {
        
        df0 = data.frame("y" = v0[[i]], "t" = v0_t, 
                         "x" = rep(0, length(v0_t)))
        df1 = data.frame("y" = v1[[i]], "t" = v1_t, 
                         "x" = rep(1, length(v1_t)))
        
        df_tot = rbind(df0, df1)
        
        init_par = beta[i,]
        
        op = optim(par=init_par, fn=min_residuals, data=df_tot, hessian = T)
        
        if(op$convergence != 0) {print("no convergence")}
        
        optim_coeff[i, ] = op$par
    }
    
    optim_coeff_split = matrix(0, nrow = 5, ncol = 3)
    for(i in 1:5) {
        
        df0 = data.frame("y" = v0_split[[i]], "t" = v0_t_split, 
                         "x" = rep(0, length(v0_t_split)))
        df1 = data.frame("y" = v1_split[[i]], "t" = v1_t_split, 
                         "x" = rep(1, length(v1_t_split)))
        
        df_tot = rbind(df0, df1)
        
        init_par = beta[i,]
        
        op = optim(par=init_par, fn=min_residuals, data=df_tot, hessian = T)
        
        if(op$convergence != 0) {print("no convergence")}
        
        optim_coeff_split[i, ] = op$par
        
        # fish_info = solve(op$hessian)
        # prop_sigma<-sqrt(diag(fish_info))
        # conf_int[i,]
    }
    
    
    
    par_est[[it]] = optim_coeff
    par_est_split[[it]] = optim_coeff_split
    
}

par_est_mat       = matrix(nrow = 100, ncol = 15)
par_est_mat_split = matrix(nrow = 100, ncol = 15)

for(i in 1:100) {
    par_est_mat[i,] = c(par_est[[i]])
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
    yVar = c(par_est_mat[,r], par_est_mat_split[,r])
    disc_type = c(rep('cond.', nrow(par_est_mat)), rep('split', nrow(par_est_mat_split)))
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