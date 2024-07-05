library(AalenJohansen)

# These initial parameters are from mcmc_runfile_deSolve.r
init_par = c(c(matrix(c( -2.348358,  0.1054161, -0.4348307,
                         -1.589180, -4.1241792,  0.7553760,
                         -1.445115, -0.0985043,  0.4890515,
                         -1.960385,  0.0327601,  0.1237441,
                         -2.533004,  0.1166776,  1.4994372), ncol=3, byrow=T)),
             c(-4.864549, -1.000305,  -2.511552,  -2.738977),
             c(-6.536839, -5.683834))

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)
beta <- matrix(init_par[par_index$beta], ncol = 3, byrow = F)

# loading the real CAV data set
load("real_cav_analysis/Data/cavData_cont.rda")
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

save(optim_coeff_split, file = paste0("real_cav_analysis/Plots/aj_out.rda"))