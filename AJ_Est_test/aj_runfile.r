# ------------------------------------------------------------------------------
# Function to calculate Aalen-Johansen estimator -------------------------------
# ------------------------------------------------------------------------------

obs_trans <- function(df) {
    count_transitions = matrix(0, nrow=4, ncol=4)
    total_trans = 0
    for(i in unique(df[,"ptnum"])){
        
        b_i_mle = c(df[df[,"ptnum"] == i, 'state'])
        
        for(t in 1:(length(b_i_mle) - 1)) {
            count_transitions[b_i_mle[t], b_i_mle[t+1]] = 
                count_transitions[b_i_mle[t], b_i_mle[t+1]] + 1
            total_trans = total_trans + 1
        }
    }
    print(count_transitions)   
}

aj_estimate <- function(s, t, data_mat, sex) {
    
    # Take only subjects with the correct covariate specification
    data_mat_sub = data_mat[data_mat$sex == sex, ]
    # data_mat_sub = data_mat
    
    # Partition the time domain
    t_unique = sort(unique(data_mat_sub$years))
    
    focus_times = t_unique[t_unique > s & t_unique <= t]
    
    A_hat = vector(mode = "list", length = length(focus_times))
    
    P_s_t = diag(4)
    
    for(t_j in focus_times) {
        print(paste0(which(focus_times == t_j), " of ", length(focus_times)))
        alpha_j = matrix(0, nrow = 4, ncol = 4)
        
        t_j_ind = which(data_mat_sub$years == t_j)
        t_j_ind_1 = t_j_ind - 1
        
        state_t_j = data_mat_sub[t_j_ind,"state"]
        state_t_j_1 = data_mat_sub[t_j_ind_1,"state"]
        
        r_j = rep(0,4)
        for(i in unique(data_mat_sub$ptnum)) {
            sub_dat = data_mat_sub[data_mat_sub$ptnum == i, ]
            t_pt = max(which(sub_dat$years < t_j))
            s_i = sub_dat[t_pt, "state"]
            r_j[s_i] = r_j[s_i] + 1
        }
        
        if(sum(r_j) != length(unique(data_mat_sub$ptnum))) print("issue with r_j")
        
        for(g in 1:4) {
            for(h in 1:4) {
                if(g != h) {
                    alpha_j[g,h] = sum(state_t_j_1 == g & state_t_j == h)
                }
            }
        }
        
        for(r in 1:nrow(alpha_j)) {
            if(r_j[r] != 0) {
                alpha_j[r, ] = alpha_j[r, ] / r_j[r]   
            }
        }
        
        diag(alpha_j) = -rowSums(alpha_j)
        
        A_hat[[which(focus_times == t_j)]] = vector(mode = 'list', length = 2)
        A_hat[[which(focus_times == t_j)]][[1]] = t_j
        A_hat[[which(focus_times == t_j)]][[2]] = alpha_j
        
        P_s_t = P_s_t %*% (diag(4) + alpha_j)
    }
    
    return_list = list(P_s_t, A_hat)
    
    return(return_list)
}

# ------------------------------------------------------------------------------
# Applying to simulated data ---------------------------------------------------
# (ignoring covariates; ignoring misclassification) ----------------------------
# ------------------------------------------------------------------------------
# it 1: no sex, slope on time is 3x
# it 2: no sex, slope on time is 1x
# it 3: sex, slope on time is 3x
# it 4: sex, slope on time is 1x

args = commandArgs(TRUE)
it = as.numeric(args[1]) 

load(paste0('DataOut/cavData', it, '.rda'))

obs_trans(cavData)

load(paste0('DataOut/trueValues.rda'))

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)

t1 = min(cavData$years)
t2 = max(cavData$years)

p_ic <- c(p1=1,p2=0,p3=0,p4=0,p5=1,p6=0,p7=0,p8=1,p9=0) # initial condition

beta <- matrix(trueValues[par_index$beta], ncol = 3, byrow = F)

beta_estimates = vector(mode = 'list', length = 2)
for(s in 0:1) {
    AJ_list = aj_estimate(t1, t2, cavData, s)
    # save(AJ_list, file = paste0("DataOut/AJ_list_sex_", s, "_", it, ".rda"))
    # load(paste0("DataOut/AJ_list_sex_", s, "_", it, ".rda"))
    
    # Store the time points and estimted transition rates and fit a regression line
    # beta_est_list[[i]][[j]] is the estimated transition rates for the i -> j transition
    beta_est_list = vector(mode = 'list', length = 4)
    for(j in 1:length(beta_est_list)) beta_est_list[[j]] = vector(mode = 'list', length = 4)
    
    for(i in 1:length(AJ_list[[2]])) {
        for(j in 1:4) {
            for(k in 1:4) {
                if(j != k) {
                    if(AJ_list[[2]][[i]][[2]][j, k] != 0) {
                        beta_est_list[[j]][[k]] = c(beta_est_list[[j]][[k]], AJ_list[[2]][[i]][[1]], AJ_list[[2]][[i]][[2]][j, k])
                    }
                }
            }
        }
    }
    
    q_comp = list()
    if(length(beta_est_list[[1]][[2]]) > 0) {
        q_comp[[1]] = t(matrix(beta_est_list[[1]][[2]], nrow = 2))   
        q_comp[[1]][,2] = log(q_comp[[1]][,2])
    } else {
        q_comp[[1]] = NULL        
    }
    if(length(beta_est_list[[1]][[4]]) > 0) {
        q_comp[[2]] = t(matrix(beta_est_list[[1]][[4]], nrow = 2))   
        q_comp[[2]][,2] = log(q_comp[[2]][,2])
    } else {
        q_comp[[2]] = NULL
    }
    if(length(beta_est_list[[2]][[3]]) > 0) {
        q_comp[[3]] = t(matrix(beta_est_list[[2]][[3]], nrow = 2))
        q_comp[[3]][,2] = log(q_comp[[3]][,2])
    } else {
        q_comp[[3]] = NULL   
    }
    if(length(beta_est_list[[2]][[4]]) > 0) {
        q_comp[[4]] = t(matrix(beta_est_list[[2]][[4]], nrow = 2))        
        q_comp[[4]][,2] = log(q_comp[[4]][,2])
    } else {
        q_comp[[4]] = NULL
    }
    if(length(beta_est_list[[3]][[4]]) > 0) {
        q_comp[[5]] = t(matrix(beta_est_list[[3]][[4]], nrow = 2))
        q_comp[[5]][,2] = log(q_comp[[5]][,2])
    } else {
        q_comp[[5]] = NULL
    }
    
    beta_estimates[[s+1]] = matrix(ncol = 2, nrow = length(q_comp))
    
    # pdf(paste0("Plots/beta_plot_sex_", s, "_", it, ".pdf"))
    # par(mfrow = c(3,2))
    for(q in 1:length(q_comp)) {
        if(!is.null(q_comp[[q]])) {
            m1 = lm(q_comp[[q]][,2] ~ q_comp[[q]][,1])
            beta_estimates[[s+1]][q, ] = c(m1$coefficients[1], m1$coefficients[2])
            
            if(s == 0) {
                title_est = paste0("beta0 = ", round(m1$coefficients[1], digits = 4),
                                   ", beta1 = ", round(m1$coefficients[2], digits = 4))
                sub_title = paste0("beta0 = ", round(beta[q,1], digits = 4),
                                   ", beta1 = ", round(beta[q,2], digits = 4))
            } else {
                title_est = paste0("beta0 = ", round(m1$coefficients[1], digits = 4),
                                   ", beta1 = ", round(m1$coefficients[2], digits = 4))
                sub_title = paste0("beta0 = ", round(beta[q,1] + beta[q,3], digits = 4),
                                   ", beta1 = ", round(beta[q,2], digits = 4))
            }
            # plot(q_comp[[q]][,1], q_comp[[q]][,2], xlab = "time",
            #      ylab = paste0("log(q",q,")"), main = title_est, sub = sub_title, col.sub = 'blue')
            # abline(m1, col = 'red')
        } else {
            # plot.new()
        }
    }
    # dev.off()
    
}

save(beta_estimates, file = paste0("DataOut/beta_estimates_", it, ".rda"))

# ------------------------------------------------------------------------------
# Comparing Numerical Integration to Product Integration -----------------------
# ------------------------------------------------------------------------------

library(deSolve, quietly=T)
Q <- function(time,sex,betaMat){

    q1  = exp( c(1,time,sex) %*% betaMat[1,] )  # Transition from state 1 to state 2.
    q2  = exp( c(1,time,sex) %*% betaMat[2,] )  # Transition from state 1 to death.
    q3  = exp( c(1,time,sex) %*% betaMat[3,] )  # Transition from state 2 to state 3.
    q4  = exp( c(1,time,sex) %*% betaMat[4,] )  # Transition from state 2 to death.
    q5  = exp( c(1,time,sex) %*% betaMat[5,] )  # Transition from state 3 to death.

    qmat = matrix(c( 0,q1, 0,q2,
                     0, 0,q3,q4,
                     0, 0, 0,q5,
                     0, 0, 0, 0),nrow=4,byrow=TRUE)
    diag(qmat) = -rowSums(qmat)

    return(qmat)
}

model_t <- function(t,p,parms) {

    betaMat <- matrix(parms$b, ncol = 3, byrow = F)

    q1  = exp( c(1,t,parms$x_ik) %*% betaMat[1,] )  # Transition from state 1 to state 2.
    q2  = exp( c(1,t,parms$x_ik) %*% betaMat[2,] )  # Transition from state 1 to death.
    q3  = exp( c(1,t,parms$x_ik) %*% betaMat[3,] )  # Transition from state 2 to state 3.
    q4  = exp( c(1,t,parms$x_ik) %*% betaMat[4,] )  # Transition from state 2 to death.
    q5  = exp( c(1,t,parms$x_ik) %*% betaMat[5,] )  # Transition from state 3 to death.

    dP = rep(1,9) # this is the vector with all differential equations

    dP[1] = p[1]*(-q1-q2)
    dP[2] = p[1]*q1 + p[2]*(-q3-q4)
    dP[3] = p[2]*q3 - p[3]*q5
    dP[4] = p[1]*q2 + p[2]*q4 + p[3]*q5
    dP[5] = p[5]*(-q3-q4)
    dP[6] = p[5]*q3 - p[6]*q5
    dP[7] = p[5]*q4 + p[6]*q5
    dP[8] = -p[8]*q5
    dP[9] = p[8]*q5

    return(list(dP))

}
                    
p_ic <- c(p1=1,p2=0,p3=0,p4=0,p5=1,p6=0,p7=0,p8=1,p9=0) # initial condition

t1 = 0
t2 = 5

load(paste0('DataOut/trueValues.rda'))
par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)
beta <- matrix(trueValues[par_index$beta], ncol = 3, byrow = F)
initProbs_temp = c( 1, exp(trueValues[par_index$pi_logit][1]), exp(trueValues[par_index$pi_logit][2]), 0)
initProbs = initProbs_temp / sum(initProbs_temp)

out1 <- deSolve::ode(p_ic, times = c(t1,t2), func = model_t,
                     parms = list(b=c(beta), x_ik = 1))

P_desolve1 <- matrix(c(out1[2,"p1"], out1[2,"p2"], out1[2,"p3"], out1[2,"p4"],
                      0, out1[2,"p5"], out1[2,"p6"], out1[2,"p7"],
                      0,  0, out1[2,"p8"], out1[2,"p9"],
                      0,  0,  0,  1), nrow = 4, byrow = T)

out0 <- deSolve::ode(p_ic, times = c(t1,t2), func = model_t,
                     parms = list(b=c(beta), x_ik = 0))

P_desolve0 <- matrix(c(out0[2,"p1"], out0[2,"p2"], out0[2,"p3"], out0[2,"p4"],
                      0, out0[2,"p5"], out0[2,"p6"], out0[2,"p7"],
                      0,  0, out0[2,"p8"], out0[2,"p9"],
                      0,  0,  0,  1), nrow = 4, byrow = T)

library(AalenJohansen)
P_prodInt_list1 <- prodint(t1, t2, 0.01, function(t){Q(t, sex = 1, betaMat = beta)})
P_prodInt_list0 <- prodint(t1, t2, 0.01, function(t){Q(t, sex = 0, betaMat = beta)})
P_prodInt1 = P_prodInt_list1[[length(P_prodInt_list1)]]
P_prodInt0 = P_prodInt_list0[[length(P_prodInt_list0)]]

print(P_desolve1); print(P_prodInt1)
print(P_desolve0); print(P_prodInt0)

# Format the simulated data into the form for AalenJohansen
par_est = vector(mode = 'list', length = 100)
for(it in 1:100) {
    print(it)
    load(paste0('DataOut/cavData', it, '.rda'))
    eid = unique(cavData$ptnum)
    cavData_aj = vector(mode = 'list', length = length(eid))
    for(i in 1:length(cavData_aj)) {
        sub_dat = cavData[cavData$ptnum == eid[i], ]
        cavData_aj[[i]] = list(times = sub_dat$years, states = sub_dat$state, X = sub_dat$sex[1])
    }
    
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
    
    if(length(v1_t) > length(v0_t)) {
        if(sum(v0_t %in% v1_t) != length(v0_t)) print(paste0(it, ": different times"))
    } else {
        if(sum(v1_t %in% v0_t) != length(v1_t)) print(paste0(it, ": different times"))
    }
    
    v_both = list()
    for(jj in 1:5) {
        df0 = data.frame("y" = v0[[jj]], "t" = v0_t, "x" = rep(0, length(v0_t)))
        df1 = data.frame("y" = v1[[jj]], "t" = v1_t, "x" = rep(1, length(v1_t)))
        df_tot = rbind(df0, df1)
        
        df_tot = df_tot[-which(df_tot$y == 0),]
        df_tot$y = log(df_tot$y)
        v_both[[jj]] = df_tot
    }
    
    coeff_sum = matrix(nrow = 5, ncol = 3)
    for(jj in 1:5) {
        if(nrow(v_both[[jj]]) > 0) {
            lm_1 = lm(y ~ t + x, data = v_both[[jj]])
            coeff_sum[jj,] = lm_1$coefficients   
        } else {
            print(paste0("No transitions for ", it, " transition ", jj))
        }
    }
    
    par_est[[it]] = coeff_sum
    
}

par_est_mat = matrix(nrow = 100, ncol = 15)
for(i in 1:100) {
    par_est_mat[i,] = c(par_est[[i]])
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
pdf(paste0('Plots/par_est_ideal.pdf'))

VP <- vector(mode="list", length = length(labels))
for(r in 1:length(labels)) {
    # Adding the boxplots
    yVar = par_est_mat[,r]
    disc_type = rep(1, nrow(par_est_mat))
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

# p1 <- unlist(lapply(fit1$p, FUN = function(L) L[2]))
# P1 <- unlist(lapply(prodint(0, 15, 0.01, function(t){Q(t, sex = 1, betaMat = beta)}),
#                     FUN = function(L) (initProbs %*% L)[2]))
# p2 <- unlist(lapply(fit0$p, FUN = function(L) L[2]))
# P2 <- unlist(lapply(prodint(0, 15, 0.01, function(t){Q(t, sex = 0, betaMat = beta)}),
#                     FUN = function(L) (initProbs %*% L)[2]))
# 
# pdf("Plots/q_AJ.pdf")
# par(mfrow = c(3, 2))
# par(mar = c(2.5, 2.5, 1.5, 1.5))
# 
# plot(v1_t, v1_12, type = "l", lty = 2, xlab = "", ylab = "", lwd = 3, 
#      main = "Hazard (1->2)", col = "red", ylim = c(min(v1_12, v0_12), max(v1_12, v0_12)))
# lines(v0_t, v0_12, lty = 2, col = "blue", lwd = 3, )
# 
# 
# plot(v1_t, v1_14, type = "l", lty = 2, xlab = "", ylab = "", lwd = 3,  
#      main = "Hazard (1->4)", col = "red", ylim = c(min(v1_14, v0_14), max(v1_14, v0_14)))
# lines(v0_t, v0_14, lty = 2, col = "blue", lwd = 3, )
# plot(v1_t, v1_23, type = "l", lty = 2, xlab = "", ylab = "", lwd = 3,  
#      main = "Hazard (2->3)", col = "red", ylim = c(min(v1_23, v0_23), max(v1_23, v0_23)))
# lines(v0_t, v0_23, lty = 2, col = "blue", lwd = 3, )
# plot(v1_t, v1_24, type = "l", lty = 2, xlab = "", ylab = "", lwd = 3,  
#      main = "Hazard (2->4)", col = "red", ylim = c(min(v1_24, v0_24), max(v1_24, v0_24)))
# lines(v0_t, v0_24, lty = 2, col = "blue", lwd = 3, )
# plot(v1_t, v1_34, type = "l", lty = 2, xlab = "", ylab = "", lwd = 3,  
#      main = "Hazard (3->4)", col = "red", ylim = c(min(v1_34, v0_34), max(v1_34, v0_34)))
# lines(v0_t, v0_34, lty = 2, col = "blue", lwd = 3, )
# 
# dev.off()

# plot(v10, p1, type = "l", lty = 2, xlab = "", ylab = "", main = "Probability", col = "red")
# lines(seq(0, 15, 0.01), P1, col = "red")
# lines(v20, p2, lty = 2, col = "blue")
# lines(seq(0, 15, 0.01), P2, col = "blue")
