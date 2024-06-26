library(AalenJohansen)
library(nhm)

args = commandArgs(TRUE)
case_num = as.numeric(args[1])

# case_num = 1 --> not exact transition time nor all transitions observed
# case_num = 2 --> exact transition times
# case_num = 3 --> exact transition times and all transitions observed

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)

# load('real_cav_analysis/Model_out/deSolve/mcmc_out_10.rda')
load('mcmc_out_10.rda')
chain = mcmc_out$chain[10000:25001, ]
ind_keep = seq(1, nrow(chain), by=10)
chain = chain[ind_keep, ]
trueValues = colMeans(chain)
trueValues[6:10] = 3 * trueValues[6:10]
trueValues[7] = trueValues[8]
trueValues[8] = trueValues[6]

beta <- matrix(trueValues[par_index$beta], ncol = 3, byrow = F)

for(it in 1:100) {
    print(it)

    # Format the simulated data into the form for AalenJohansen
    load(paste0("DataOut/cavData_case", case_num, "_it", it, ".rda"))
    
    eid = unique(cavData$ptnum)
    
    cavData_aj = vector(mode = 'list', length = length(eid))
    for(i in 1:length(cavData_aj)) {
        sub_dat = cavData[cavData$ptnum == eid[i], ]
        cavData_aj[[i]] = list(times = sub_dat$years, states = sub_dat$state, X = sub_dat$sex[1])
    }
    
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
    
    #  -------------------------------------------------------------------------
    # Using optim to get parameter estimates -----------------------------------
    #  -------------------------------------------------------------------------
    min_residuals <- function(data, par) {
        with(data, sum(abs( (1/par[2]) * (exp(par[1] + par[2]*t + par[3]*x) - exp(par[1] + par[3]*x)) - y)))
    }
    
    optim_coeff_split = matrix(0, nrow = 5, ncol = 3)
    for(i in 1:5) {
        
        df0 = data.frame("y" = v0_split[[i]], "t" = v0_t_split, 
                         "x" = rep(0, length(v0_t_split)))
        df1 = data.frame("y" = v1_split[[i]], "t" = v1_t_split, 
                         "x" = rep(1, length(v1_t_split)))
        
        df_tot = rbind(df0, df1)
        
        if(sum(!is.finite(df_tot$y)) > 0) {
            print(paste0(it, ', ', i, ', has inf'))
            df_tot = df_tot[is.finite(df_tot$y), ]
        }
        
        init_par = beta[i,]
        
        op = optim(par=init_par, fn=min_residuals, data=df_tot, control = list(maxit = 1e5))
        
        if(op$convergence != 0) {print("no convergence")}
        
        optim_coeff_split[i, ] = op$par
    }
    
    par_est_list = list(c(optim_coeff_split))
    
    save(par_est_list, file = paste0('Model_out/par_est_list_case', case_num, "_it", it, '.rda'))
}
