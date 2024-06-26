library(AalenJohansen)

args = commandArgs(TRUE)
exact_time = as.logical(as.numeric(args[1]))

for(it in 1:100) {
    print(it)
    
    par_index = list( beta=1:18)
    
    trueValues=c(matrix(c(-3,  0.50, 0.75,
                          -3,  0.50, 0.75,
                          -3,  0.50, 0.75,
                          -3,  0.50, 0.75,
                          -3,  0.50, 0.75,
                          -3,  0.50, 0.75), ncol = 3, byrow = T))
    
    beta <- matrix(trueValues[par_index$beta], ncol = 3, byrow = F)
    
    # Format the simulated data into the form for AalenJohansen
    if(exact_time) {
        load(paste0('DataOut/exactTime/cavData', it, '.rda'))
    } else {
        load(paste0('DataOut/interTime/cavData', it, '.rda'))
    }
    
    # Sample size of 2000
    cavData = cavData[cavData$ptnum %in% c(1:2000), ]
    
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
    v1_split[[2]] <- unlist(lapply(fit1_split$Lambda, FUN = function(L) L[1,3]))
    v1_split[[3]] <- unlist(lapply(fit1_split$Lambda, FUN = function(L) L[2,1]))
    v1_split[[4]] <- unlist(lapply(fit1_split$Lambda, FUN = function(L) L[2,3]))
    v1_split[[5]] <- unlist(lapply(fit1_split$Lambda, FUN = function(L) L[3,1]))
    v1_split[[6]] <- unlist(lapply(fit1_split$Lambda, FUN = function(L) L[3,2]))
    
    v0_split = list()
    v0_split[[1]] <- unlist(lapply(fit0_split$Lambda, FUN = function(L) L[1,2]))
    v0_split[[2]] <- unlist(lapply(fit0_split$Lambda, FUN = function(L) L[1,3]))
    v0_split[[3]] <- unlist(lapply(fit0_split$Lambda, FUN = function(L) L[2,1]))
    v0_split[[4]] <- unlist(lapply(fit0_split$Lambda, FUN = function(L) L[2,3]))
    v0_split[[5]] <- unlist(lapply(fit0_split$Lambda, FUN = function(L) L[3,1]))
    v0_split[[6]] <- unlist(lapply(fit0_split$Lambda, FUN = function(L) L[3,2]))
    
    v1_t_split <- fit1_split$t
    v0_t_split <- fit0_split$t
    
    #  -------------------------------------------------------------------------
    # Using optim to get parameter estimates -----------------------------------
    #  -------------------------------------------------------------------------
    min_residuals <- function(data, par) {
        with(data, sum(( (1/par[2]) * (exp(par[1] + par[2]*t + par[3]*x) - exp(par[1] + par[3]*x)) - y)^2))
    }
    
    optim_coeff_split = matrix(0, nrow = 6, ncol = 3)
    for(i in 1:nrow(optim_coeff_split)) {
        
        df0 = data.frame("y" = v0_split[[i]], "t" = v0_t_split, 
                         "x" = rep(0, length(v0_t_split)))
        df1 = data.frame("y" = v1_split[[i]], "t" = v1_t_split, 
                         "x" = rep(1, length(v1_t_split)))
        
        df_tot = rbind(df0, df1)
        
        if(sum(!is.finite(df_tot$y)) > 0) {
            print(paste0(it, ', ', i, ', has inf'))
            df_tot = df_tot[is.finite(df_tot$y), ]
        }
        
        init_par = beta[i,] + runif(n = 3, min = -0.25, max = 0.25)
        # init_par = beta[i,]
        
        op = optim(par=init_par, fn=min_residuals, data=df_tot, control = list(maxit = 1e5)) #reltol = 1e-5
        
        if(op$convergence != 0) {print("no convergence")}
        
        optim_coeff_split[i, ] = op$par
    }
    
    par_est_list = list(c(optim_coeff_split))
    
    if(exact_time) {
        save(par_est_list, file = paste0('Model_out/exactTime/par_est_list_', it, '.rda'))
    } else {
        save(par_est_list, file = paste0('Model_out/interTime/par_est_list_', it, '.rda'))
    }    
}