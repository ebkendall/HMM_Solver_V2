# ------------------------------------------------------------------------------
# Using the AalenJohansen estimator package ------------------------------------
# ------------------------------------------------------------------------------
library(AalenJohansen)

t1 = 0
t2 = 5

load(paste0('DataOut/trueValues_4.rda'))
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

P_prodInt_list1 <- prodint(t1, t2, 0.01, function(t){Q(t, sex = 1, betaMat = beta)})
P_prodInt_list0 <- prodint(t1, t2, 0.01, function(t){Q(t, sex = 0, betaMat = beta)})
P_prodInt1 = P_prodInt_list1[[length(P_prodInt_list1)]]
P_prodInt0 = P_prodInt_list0[[length(P_prodInt_list0)]]

print(P_desolve1); print(P_prodInt1)
print(P_desolve0); print(P_prodInt0)

# Format the simulated data into the form for AalenJohansen
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