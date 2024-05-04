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
    
    # Partition the time domain
    t_unique = sort(unique(data_mat$years))
    
    # Take only subjects with the correct covariate specification
    data_mat_sub = data_mat[data_mat$sex == sex, ]
    
    focus_times = t_unique[t_unique > s & t_unique <= t]
    
    A_hat = vector(mode = "list", length = length(focus_times))
    
    P_s_t = diag(4)
    
    for(t_j in focus_times) {
        print(paste0(which(focus_times == t_j), " of ", length(focus_times)))
        alpha_j = matrix(0, nrow = 4, ncol = 4)
        
        if(t_j %in% data_mat_sub$years) {
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
        } else {
            A_hat[[which(focus_times == t_j)]] = vector(mode = 'list', length = 2)
            A_hat[[which(focus_times == t_j)]][[1]] = t_j
            A_hat[[which(focus_times == t_j)]][[2]] = alpha_j
        }
        
        
        P_s_t = P_s_t %*% (diag(4) + alpha_j)
    }
    
    return_list = list(P_s_t, A_hat)
    
    return(return_list)
}

aj_estimate_new <- function(data_mat, sex) {
    
    data_mat_sub = data_mat[data_mat[,"sex"] == sex, ]
    
    s = 0
    t = max(data_mat_sub[,"years"])
    
    EIDs = unique(data_mat_sub[,"ptnum"])
    
    t_j = NULL      # times when transitions between states are recorded
    t_j_ind = NULL  # row index of data_mat_sub where each t_j lives
    
    for(i in 1:length(EIDs)) {
        t_j_ind_i = which(data_mat_sub[,"ptnum"] == EIDs[i])
        df_i = data_mat_sub[t_j_ind_i, ]
        
        s_change = diff(df_i[,"state"])
        
        s_change_pos = which(s_change != 0) + 1
        
        # Check to see if we observed any state changes
        if(length(s_change_pos) > 0) {
            t_j = c(t_j, data_mat_sub[t_j_ind_i[s_change_pos],"years"])
            t_j_ind = c(t_j_ind, t_j_ind_i[s_change_pos])
        }
    }
    
    t_mat = cbind(t_j, t_j_ind)
    
    unique_times = sort(unique(t_j))
    
    alpha_hat_list = vector(mode = 'list', length = length(unique_times))
    
    track_total_transitions = 0
    
    for(j in 1:length(unique_times)) {
        
        print(paste0(j, " of ", length(unique_times)))
        
        t_mat_sub = t_mat[t_mat[,1] == unique_times[j], , drop = F]
        
        d_mat_j = matrix(0, nrow = 4, ncol = 4)
        r_j = rep(0, 4)
        
        # Calculate d_ghj
        for(k in 1:nrow(t_mat_sub)) {
            s_to_ind = t_mat_sub[k, 2]
            s_from_ind = s_to_ind - 1
            
            s_to = data_mat_sub[s_to_ind, "state"]
            s_from = data_mat_sub[s_from_ind, "state"]
            
            d_mat_j[s_from, s_to] = d_mat_j[s_from, s_to] + 1
            track_total_transitions = track_total_transitions + 1
        }
        
        # Calculate r_gj
        data_mat_prior_t_j = data_mat_sub[data_mat_sub[,"years"] < unique_times[j], ,drop = F]
        sub_prio_t_j = unique(data_mat_prior_t_j[,"ptnum"])
        for(u in 1:length(sub_prio_t_j)) {
            df_t_j = data_mat_prior_t_j[data_mat_prior_t_j[,"ptnum"] == sub_prio_t_j[u], , drop = F]
            curr_state = df_t_j[nrow(df_t_j), "state"]
            r_j[curr_state] = r_j[curr_state] + 1
        }
        
        if(sum(r_j) != length(EIDs)) {print("something went wrong (1)")}
        
        alpha_hat_list[[j]] = list()
        alpha_hat_list[[j]][[1]] = unique_times[j]
        alpha_hat_list[[j]][[2]] = d_mat_j
        alpha_hat_list[[j]][[3]] = r_j
    }
    
    if(track_total_transitions != length(t_j)) {
        print("something went wrong (2)")
        print(paste0("total transitions = ", track_total_transitions))
        print(paste0("total t_j's = ", length(t_j)))
    }
    
    return(alpha_hat_list)
}

# args = commandArgs(TRUE)
# it = as.numeric(args[1]) 
it = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

load(paste0('DataOut/cavData', it, '.rda'))

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)

t1 = min(cavData$years)
t2 = max(cavData$years)

beta_estimates = vector(mode = 'list', length = 2)
for(s in 0:1) {
    AJ_list = aj_estimate(t1, t2, cavData, s)
    
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
    
    beta_estimates[[s+1]] = q_comp
    
}

# Restructure to join the two groups
v_both = list()
for(jj in 1:5) {
    df0 = data.frame("y" = beta_estimates[[1]][[jj]][,2], 
                     "t" = beta_estimates[[1]][[jj]][,1], 
                     "x" = rep(0, nrow(beta_estimates[[1]][[jj]])))
    df1 = data.frame("y" = beta_estimates[[2]][[jj]][,2], 
                     "t" = beta_estimates[[2]][[jj]][,1], 
                     "x" = rep(1, nrow(beta_estimates[[2]][[jj]])))
    
    df_tot = rbind(df0, df1)
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

beta_estimates[[3]] = coeff_sum

save(beta_estimates, file = paste0("DataOut/beta_estimates_", it, ".rda"))

#------------------------------------------------------------------------------
# Testing new function --------------------------------------------------------
#------------------------------------------------------------------------------
alpha_hat_list1 = aj_estimate_new(cavData, 1)
alpha_hat_list0 = aj_estimate_new(cavData, 0)

# sex = 1
q1_hat = vector(mode = 'list', length = 5)
for(j in 1:length(alpha_hat_list1)) {
    
    if(alpha_hat_list1[[j]][[2]][1,2] != 0) {
        t = alpha_hat_list1[[j]][[1]]
        q = alpha_hat_list1[[j]][[2]][1,2] / alpha_hat_list1[[j]][[3]][1]
        
        q1_hat[[1]] = rbind(q1_hat[[1]], c(t,q))
    }
    
    if(alpha_hat_list1[[j]][[2]][1,4] != 0) {
        t = alpha_hat_list1[[j]][[1]]
        q = alpha_hat_list1[[j]][[2]][1,4] / alpha_hat_list1[[j]][[3]][1]
        
        q1_hat[[2]] = rbind(q1_hat[[2]], c(t,q))
    }
    
    if(alpha_hat_list1[[j]][[2]][2,3] != 0) {
        t = alpha_hat_list1[[j]][[1]]
        q = alpha_hat_list1[[j]][[2]][2,3] / alpha_hat_list1[[j]][[3]][2]
        
        q1_hat[[3]] = rbind(q1_hat[[3]], c(t,q))
    }

    if(alpha_hat_list1[[j]][[2]][2,4] != 0) {
        t = alpha_hat_list1[[j]][[1]]
        q = alpha_hat_list1[[j]][[2]][2,4] / alpha_hat_list1[[j]][[3]][2]
        
        q1_hat[[4]] = rbind(q1_hat[[4]], c(t,q))
    }
    
    if(alpha_hat_list1[[j]][[2]][3,4] != 0) {
        t = alpha_hat_list1[[j]][[1]]
        q = alpha_hat_list1[[j]][[2]][3,4] / alpha_hat_list1[[j]][[3]][3]
        
        q1_hat[[5]] = rbind(q1_hat[[5]], c(t,q))
    }
}

# sex = 0
q0_hat = vector(mode = 'list', length = 5)
for(j in 1:length(alpha_hat_list0)) {
    
    if(alpha_hat_list0[[j]][[2]][1,2] != 0) {
        t = alpha_hat_list0[[j]][[1]]
        q = alpha_hat_list0[[j]][[2]][1,2] / alpha_hat_list0[[j]][[3]][1]
        
        q0_hat[[1]] = rbind(q0_hat[[1]], c(t,q))
    }
    
    if(alpha_hat_list0[[j]][[2]][1,4] != 0) {
        t = alpha_hat_list0[[j]][[1]]
        q = alpha_hat_list0[[j]][[2]][1,4] / alpha_hat_list0[[j]][[3]][1]
        
        q0_hat[[2]] = rbind(q0_hat[[2]], c(t,q))
    }
    
    if(alpha_hat_list0[[j]][[2]][2,3] != 0) {
        t = alpha_hat_list0[[j]][[1]]
        q = alpha_hat_list0[[j]][[2]][2,3] / alpha_hat_list0[[j]][[3]][2]
        
        q0_hat[[3]] = rbind(q0_hat[[3]], c(t,q))
    }
    
    if(alpha_hat_list0[[j]][[2]][2,4] != 0) {
        t = alpha_hat_list0[[j]][[1]]
        q = alpha_hat_list0[[j]][[2]][2,4] / alpha_hat_list0[[j]][[3]][2]
        
        q0_hat[[4]] = rbind(q0_hat[[4]], c(t,q))
    }
    
    if(alpha_hat_list0[[j]][[2]][3,4] != 0) {
        t = alpha_hat_list0[[j]][[1]]
        q = alpha_hat_list0[[j]][[2]][3,4] / alpha_hat_list0[[j]][[3]][3]
        
        q0_hat[[5]] = rbind(q0_hat[[5]], c(t,q))
    }
}


n_trans1 = obs_trans(cavData[cavData$sex == 1, ])
n_trans0 = obs_trans(cavData[cavData$sex == 0, ])

#------------------------------------------------------------------------------
# Plotting the results --------------------------------------------------------
#------------------------------------------------------------------------------
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

labels = c("S1->S2 (sex = 1)", "S1->S4 (sex = 1)", "S2->S3 (sex = 1)",
           "S2->S4 (sex = 1)", "S3->S4 (sex = 1)",
           "S1->S2 (sex = 0)", "S1->S4 (sex = 0)", "S2->S3 (sex = 0)",
           "S2->S4 (sex = 0)", "S3->S4 (sex = 0)")

load('DataOut/trueValues.rda')
par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)
beta <- matrix(trueValues[par_index$beta], ncol = 3, byrow = F)

r1 = list()
r0 = list()
rate_mat_ind = c(5, 13, 10, 14, 15)
for(i in 1:5) {
    r1_times = q1_hat[[i]][,1]
    r0_times = q0_hat[[i]][,1]
    
    r1[[i]] = matrix(nrow = nrow(q1_hat[[i]]), ncol = 2); r1[[i]][,1] = r1_times
    r0[[i]] = matrix(nrow = nrow(q0_hat[[i]]), ncol = 2); r0[[i]][,1] = r0_times
    
    for(t in 1:length(r1_times)) {
        r1[[i]][t,2] = c(Q(r1_times[t], 1, beta))[rate_mat_ind[i]]
    }
    for(t in 1:length(r0_times)) {
        r0[[i]][t,2] = c(Q(r0_times[t], 0, beta))[rate_mat_ind[i]]
    }
}

pdf("Plots/rate_plot_EK.pdf")
par(mfrow=c(3,2))
for(i in 1:5) {
    plot(q1_hat[[i]][,1], log(q1_hat[[i]][,2]), type = "l", lwd = 1, 
         xlab = "", ylab = "", main = labels[i], col = "red")
    lines(r1[[i]][,1], log(r1[[i]][,2]), lty = 1, col = "black")
    if(i == 1) {
        legend(x = "topleft",          # Position
               legend = c("AJ", "true"),  # Legend texts
               lty = c(1,1),           # Line types
               col = c("red", "black"),           # Line colors
               lwd = c(1, 1))
    }
    
    plot(q0_hat[[i]][,1], log(q0_hat[[i]][,2]), type = "l", lwd = 1, 
         xlab = "", ylab = "", main = labels[i+5], col = "red")
    lines(r0[[i]][,1], log(r0[[i]][,2]), lty = 1, col = "black")
    
}
dev.off()

# par_est_mat = matrix(nrow = 100, ncol = 15)
# for(i in 1:100) {
#     load(paste0("DataOut/beta_estimates_", i, ".rda"))
#     par_est_mat[i,] = c(beta_estimates[[3]])
# }
# 
# # Plot and save the mcmc trace plots and histograms.
# library(tidyverse)
# library(gridExtra)
# library(latex2exp)
# pdf(paste0('Plots/par_est_EK.pdf'))
# 
# VP <- vector(mode="list", length = length(labels))
# for(r in 1:length(labels)) {
#     # Adding the boxplots
#     yVar = par_est_mat[,r]
#     disc_type = rep(1, nrow(par_est_mat))
#     x_label = paste0("Parameter Value: ", round(trueValues[r], 3))
#     truth_par = trueValues[r]
#     
#     plot_df = data.frame(yVar = yVar, disc_type = disc_type)
#     VP[[r]] = ggplot(plot_df, aes(x=disc_type, y = yVar)) +
#         geom_violin(trim=FALSE) +
#         geom_boxplot(width=0.1) +
#         ggtitle(labels[r]) +
#         ylab(" ") +
#         xlab(x_label) +
#         geom_hline(yintercept=truth_par, linetype="dashed", color = "red") +
#         theme(text = element_text(size = 7))
# }
# 
# grid.arrange(VP[[1]], VP[[2]], VP[[3]], VP[[4]], VP[[5]], ncol=2, nrow =3)
# grid.arrange(VP[[6]], VP[[7]], VP[[8]], VP[[9]], VP[[10]], ncol=2, nrow =3)
# grid.arrange(VP[[11]], VP[[12]], VP[[13]], VP[[14]], VP[[15]], ncol=2, nrow =3)
# 
# dev.off()

