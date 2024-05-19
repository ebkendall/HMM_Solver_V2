library(AalenJohansen)
library(nhm)

# args = commandArgs(TRUE)
# it = as.numeric(args[1])
# exact_time = as.logical(as.numeric(args[2]))
# it = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
exact_time = T

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)

trueValues= c(-2.31617310,  -1.28756312,  -1.10116400,  -2.52367543,  -2.10384797,
              0.27050001, -11.65470594,  -0.49306415,   0.28862090,   0.22731278,
              -0.39079609,  -0.05894252,  -0.32509646,   0.48631653,   0.99565810,
              -5.28923943,  -0.90870027,  -2.40751854,  -2.44696544,  -6.52252202,
              -6.24090500)

beta <- matrix(trueValues[par_index$beta], ncol = 3, byrow = F)

initProbs = c(1,0,0,0)

for(it in 1:100) {
    print(it)
    # Format the simulated data into the form for AalenJohansen
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
    
    #  -------------------------------------------------------------------------
    # Using nhm to compare parameter estimates ---------------------------------
    #  -------------------------------------------------------------------------
    # nhm_data = cavData[, c("state", "years", "ptnum", "sex")]
    # colnames(nhm_data) = c("state", "time", "id", "cov1")
    # nhm_data = as.data.frame(nhm_data)
    # nhm_data$state = as.numeric(nhm_data$state)
    # 
    # trans <- rbind(c(0, 1, 0, 2),
    #                c(0, 0, 3, 4),
    #                c(0, 0, 0, 5),
    #                rep(0,4))
    # 
    # nonh <- rbind(c(0, 1, 0, 2),
    #               c(0, 0, 3, 4),
    #               c(0, 0, 0, 5),
    #               rep(0,4))
    # 
    # covm <- list("cov1" = rbind(c(0, 1, 0, 2),
    #                             c(0, 0, 3, 4),
    #                             c(0, 0, 0, 5),
    #                             rep(0,4)))
    # 
    # gomp_model <- model.nhm(state~time, data=nhm_data, subject = id, 
    #                         covariates="cov1",type="gompertz",trans=trans,
    #                         nonh=nonh,covm=covm,death = T, death.states = 4, centre_time = 5)
    # init_par <- trueValues[par_index$beta]
    # init_par[1:5] <- init_par[1:5] + 5*init_par[6:10]
    # 
    # split_max = floor(max(cavData$years))
    # 
    # gomp_fit  <- nhm(gomp_model, initial=init_par, 
    #                  control=nhm.control(obsinfo = FALSE, splits=c(seq(0.2,split_max,by=0.2))))
    # nhm_par_est = gomp_fit$par
    # nhm_par_est[1:5] = nhm_par_est[1:5] - 5*nhm_par_est[6:10]
    
    # par_est_list = list(c(optim_coeff_split), nhm_par_est)
    par_est_list = list(c(optim_coeff_split))
    
    if(exact_time) {
        save(par_est_list, file = paste0('Model_out/exactTime/par_est_list_', it, '.rda'))
    } else {
        save(par_est_list, file = paste0('Model_out/interTime/par_est_list_', it, '.rda'))
    }   
}
