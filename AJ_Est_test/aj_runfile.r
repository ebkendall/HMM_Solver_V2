library(AalenJohansen)
library(nhm)

args = commandArgs(TRUE)
it = as.numeric(args[1])
exact_time = as.logical(as.numeric(args[2]))
# it = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# exact_time = T

print(it)

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)

trueValues= c(-2.31617310,  -1.28756312,  -1.10116400,  -2.52367543,  -2.10384797,
              0.27050001, -11.65470594,  -0.49306415,   0.28862090,   0.22731278,
              -0.39079609,  -0.05894252,  -0.32509646,   0.48631653,   0.99565810,
              -5.28923943,  -0.90870027,  -2.40751854,  -2.44696544,  -6.52252202,
              -6.24090500)

beta <- matrix(trueValues[par_index$beta], ncol = 3, byrow = F)

initProbs = c(1,0,0,0)

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
nhm_data = cavData[, c("state", "years", "ptnum", "sex")]
colnames(nhm_data) = c("state", "time", "id", "cov1")
nhm_data = as.data.frame(nhm_data)
nhm_data$state = as.numeric(nhm_data$state)

trans <- rbind(c(0, 1, 0, 2),
               c(0, 0, 3, 4),
               c(0, 0, 0, 5),
               rep(0,4))

nonh <- rbind(c(0, 1, 0, 2),
              c(0, 0, 3, 4),
              c(0, 0, 0, 5),
              rep(0,4))

covm <- list("cov1" = rbind(c(0, 1, 0, 2),
                            c(0, 0, 3, 4),
                            c(0, 0, 0, 5),
                            rep(0,4)))

gomp_model <- model.nhm(state~time, data=nhm_data, subject = id, 
                        covariates="cov1",type="gompertz",trans=trans,
                        nonh=nonh,covm=covm,death = T, death.states = 4, centre_time = 5)
init_par <- trueValues[par_index$beta]
init_par[1:5] <- init_par[1:5] + 5*init_par[6:10]

split_max = floor(max(cavData$years))

gomp_fit  <- nhm(gomp_model, initial=init_par, 
                 control=nhm.control(obsinfo = FALSE, splits=c(seq(0.2,split_max,by=0.2))))
nhm_par_est = gomp_fit$par
nhm_par_est[1:5] = nhm_par_est[1:5] - 5*nhm_par_est[6:10]

par_est_list = list(c(optim_coeff_split), nhm_par_est)

if(exact_time) {
    save(par_est_list, file = paste0('Model_out/exactTime/par_est_list_', it, '.rda'))
} else {
    save(par_est_list, file = paste0('Model_out/interTime/par_est_list_', it, '.rda'))
}

# # Create a coarsened version of the dataset for numerical stability --------
# nhm_data2 <- nhm_data
# nhm_data2$time <- floor(nhm_data2$time*1000)/1000 #Coarsen time
# timegap<-c(0,diff(nhm_data2$time))
# nhm_data2$time[timegap < 1e-10 & nhm_data2$time!=0] <- nhm_data2$time[timegap<1e-10 & nhm_data2$time!=0]+0.001
# nhm_data2$time <- floor(nhm_data2$time*1000)/1000 
# gomp_model_bespoke <- model.nhm(state~time, data=nhm_data2, subject = id, 
#                                 trans=trans, covariates="cov1",
#                                 type="bespoke", intens=fourstate_illness_death,
#                                 death = T, death.states = 4)
# init_par3 <- trueValues[par_index$beta]
# init_par3[1:5] <- init_par3[1:5] + 5*init_par3[6:10]
# 
# ll<-rep(-Inf,21) 
# ul<-rep(5,21)
# 
# gomp_fit  <- nhm(gomp_model_bespoke, initial=init_par3, 
#                  control=nhm.control(obsinfo=FALSE,splits=c(seq(0.5,12,by=0.5)),
#                                      constrained=TRUE,lower_lim=ll,upper_lim=ul,
#                                      nlminb_control = list(trace=1)))
# nhm_par_est = gomp_fit$par
# nhm_par_est[1:5] = nhm_par_est[1:5] - 5*nhm_par_est[6:10]
# 
# # ------------------------------------------------------------------------------
# # nhm: Bespoke version of gompertz model, centered at 5 years: -----------------
# fourstate_illness_death<-function(t,z,x) {
#     tc <- t-5
#     q12<-exp(x[1])
#     q14<-exp(x[2])
#     q23<-exp(x[3])
#     q24<-exp(x[4])
#     q34<-exp(x[5])
#     i12<-q12*exp(x[6]*tc + x[11]*z[1])
#     i14<-q14*exp(x[7]*tc + x[12]*z[1])
#     i23<-q23*exp(x[8]*tc + x[13]*z[1])
#     i24<-q24*exp(x[9]*tc + x[14]*z[1])
#     i34<-q34*exp(x[10]*tc + x[15]*z[1])
#     q<-rbind(c(0,i12,0,i14),c(0,0,i23,i24),c(0,0,0,i34),c(0,0,0,0))
#     diag(q) <- c(-i12-i14,-i23-i24,-i34,0)
#     der<-array(0,c(4,4,15))
#     der[1,1,1]<--i12
#     der[1,2,1]<-i12
#     der[1,1,2]<--i14
#     der[1,4,2]<-i14
#     der[2,2,3]<--i23
#     der[2,3,3]<-i23
#     der[2,2,4]<--i24
#     der[2,4,4]<-i24
#     der[3,3,5]<--i34
#     der[3,4,5]<-i34
#     der[,,6:10]<-tc*der[,,1:5]
#     der[1,1,11]<--i12*z[1]
#     der[1,2,11]<-i12*z[1]
#     der[1,1,12]<--i14*z[1]
#     der[1,4,12]<-i14*z[1]
#     der[2,2,13]<--i23*z[1]
#     der[2,3,13]<-i23*z[1]
#     der[2,2,14]<--i24*z[1]
#     der[2,4,14]<-i24*z[1]
#     der[3,3,15]<--i34*z[1]
#     der[3,4,15]<-i34*z[1]
#     Q<-list(q=q,qp=der)
#     return(Q)
# }
# attr(fourstate_illness_death,"npar")<-15
# attr(fourstate_illness_death,"parnames")<-c("1->2 base:","1->4 base:","2->3 base:","2->4 base:","3->4 base:","1->2 NH","1->4 NH",
#                                             "2->3 NH:","2->4 NH:","3->4 NH:","1->2 Cov1:","1->4 Cov1:","2->3 Cov1:","2->4 Cov1:","3->4 Cov1:")
# attr(fourstate_illness_death,"parclass")<-rep(c("Trans","Nonhom","Cov"),c(5,5,5))
# # ------------------------------------------------------------------------------

