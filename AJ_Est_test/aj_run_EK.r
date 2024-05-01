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
