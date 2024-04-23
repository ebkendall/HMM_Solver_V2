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

load(paste0('DataOut/trueValues_', it, '.rda'))

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)

t1 = min(cavData$years)
t2 = max(cavData$years)

p_ic <- c(p1=1,p2=0,p3=0,p4=0,p5=1,p6=0,p7=0,p8=1,p9=0) # initial condition

beta <- matrix(trueValues[par_index$beta], ncol = 3, byrow = F)

for(s in 0:1) {
    AJ_list = aj_estimate(t1, t2, cavData, s)
    save(AJ_list, file = paste0("DataOut/AJ_list_sex_", s, "_", it, ".rda"))   
    
    # Store the time points and estimted transition rates and fit a regression line
    # beta_est_list[[i]][[j]] is the estimated transition rates for the i -> j transition
    beta_est_list = vector(mode = 'list', length = 4)
    for(j in 1:length(beta_est_list)) beta_est_list[[j]] = vector(mode = 'list', length = 4)
    
    for(i in 1:length(AJ_list[[2]])) {
        print(i)
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
    q_comp[[1]] = t(matrix(beta_est_list[[1]][[2]], nrow = 2))
    q_comp[[2]] = t(matrix(beta_est_list[[1]][[4]], nrow = 2))
    q_comp[[3]] = t(matrix(beta_est_list[[2]][[3]], nrow = 2))
    q_comp[[4]] = t(matrix(beta_est_list[[2]][[4]], nrow = 2))
    q_comp[[5]] = t(matrix(beta_est_list[[3]][[4]], nrow = 2))
    
    for(q in 1:length(q_comp)) {
        q_comp[[q]][,2] = log(q_comp[[q]][,2])
    }
    
    pdf(paste0("beta_plot_sex_", s, "_", it, ".pdf"))
    par(mfrow = c(3,2))
    for(q in 1:length(q_comp)) {
        m1 = lm(q_comp[[q]][,2] ~ q_comp[[q]][,1])
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
        plot(q_comp[[q]][,1], q_comp[[q]][,2], xlab = "time", 
             ylab = paste0("log(q",q,")"), main = title_est, sub = sub_title, col.sub = 'blue')
        abline(m1, col = 'red')   
    }
    dev.off()
    
}