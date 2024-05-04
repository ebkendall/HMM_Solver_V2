# ------------------------------------------------------------------------------
# Using the AalenJohansen estimator package ------------------------------------
# ------------------------------------------------------------------------------
library(AalenJohansen)
# ------------------------------------------------------------------------------
# Numerically solving Kolmogorov forward equations & product integral ----------
# ------------------------------------------------------------------------------
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

t1 = 0
t2 = 2

p_ic <- c(p1=1,p2=0,p3=0,p4=0,p5=1,p6=0,p7=0,p8=1,p9=0) # initial condition

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

step_size = (t2-t1)/100

P_prodInt_list1 <- prodint(t1, t2, step_size, function(t){Q(t, sex = 1, betaMat = beta)})
P_prodInt_list0 <- prodint(t1, t2, step_size, function(t){Q(t, sex = 0, betaMat = beta)})
P_prodInt1 = P_prodInt_list1[[length(P_prodInt_list1)]]
P_prodInt0 = P_prodInt_list0[[length(P_prodInt_list0)]]

print(P_desolve1); print(P_prodInt1)
print(P_desolve0); print(P_prodInt0)

# ------------------------------------------------------------------------------
# Fitting the AalenJohansen function -------------------------------------------
# ------------------------------------------------------------------------------

# Format the simulated data into the form for AalenJohansen
it = 1
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

p1 = list()
p1[[1]] <- unlist(lapply(fit1$p, FUN = function(L) L[1]))
p1[[2]] <- unlist(lapply(fit1$p, FUN = function(L) L[2]))
p1[[3]] <- unlist(lapply(fit1$p, FUN = function(L) L[3]))
p1[[4]] <- unlist(lapply(fit1$p, FUN = function(L) L[4]))

p0 = list()
p0[[1]] <- unlist(lapply(fit0$p, FUN = function(L) L[1]))
p0[[2]] <- unlist(lapply(fit0$p, FUN = function(L) L[2]))
p0[[3]] <- unlist(lapply(fit0$p, FUN = function(L) L[3]))
p0[[4]] <- unlist(lapply(fit0$p, FUN = function(L) L[4]))

v1_t <- fit1$t
v0_t <- fit0$t

P1 = list()
P1[[1]] <- unlist(lapply(prodint(0, max(v1_t), 0.01, 
                                 function(t){Q(t, sex = 1, betaMat = beta)}),
                        FUN = function(L) (initProbs %*% L)[1]))
P1[[2]] <- unlist(lapply(prodint(0, max(v1_t), 0.01, 
                                 function(t){Q(t, sex = 1, betaMat = beta)}),
                        FUN = function(L) (initProbs %*% L)[2]))
P1[[3]] <- unlist(lapply(prodint(0, max(v1_t), 0.01, 
                                 function(t){Q(t, sex = 1, betaMat = beta)}),
                        FUN = function(L) (initProbs %*% L)[3]))
P1[[4]] <- unlist(lapply(prodint(0, max(v1_t), 0.01, 
                                 function(t){Q(t, sex = 1, betaMat = beta)}),
                        FUN = function(L) (initProbs %*% L)[4]))

P0 = list()
P0[[1]] <- unlist(lapply(prodint(0, max(v0_t), 0.01, 
                                 function(t){Q(t, sex = 0, betaMat = beta)}),
                         FUN = function(L) (initProbs %*% L)[1]))
P0[[2]] <- unlist(lapply(prodint(0, max(v0_t), 0.01, 
                                 function(t){Q(t, sex = 0, betaMat = beta)}),
                         FUN = function(L) (initProbs %*% L)[2]))
P0[[3]] <- unlist(lapply(prodint(0, max(v0_t), 0.01, 
                                 function(t){Q(t, sex = 0, betaMat = beta)}),
                         FUN = function(L) (initProbs %*% L)[3]))
P0[[4]] <- unlist(lapply(prodint(0, max(v0_t), 0.01, 
                                 function(t){Q(t, sex = 0, betaMat = beta)}),
                         FUN = function(L) (initProbs %*% L)[4]))

# Gauge model fit with integrated intensity function
integrated_Q <- function(time, sex, beta) {
    p1 = exp(beta[1] + beta[2] * time + beta[3] * sex)
    p2 = exp(beta[1] + beta[3] * sex)
    
    i = (1/beta[2]) * (p1 - p2)
    return(i)
}

rate_mat1 = int_rate_mat1 = matrix(0, nrow = length(v1_t), ncol = 5)
rate_mat0 = int_rate_mat0 = matrix(0, nrow = length(v0_t), ncol = 5)
for(t in 1:length(v1_t)) {
    rate_mat1[t,] = c(Q(v1_t[t], 1, beta))[c(5, 13, 10, 14, 15)]
    
    int_rate_mat1[t,1] = integrated_Q(v1_t[t], 1, beta[1,])
    int_rate_mat1[t,2] = integrated_Q(v1_t[t], 1, beta[2,])
    int_rate_mat1[t,3] = integrated_Q(v1_t[t], 1, beta[3,])
    int_rate_mat1[t,4] = integrated_Q(v1_t[t], 1, beta[4,])
    int_rate_mat1[t,5] = integrated_Q(v1_t[t], 1, beta[5,])
}
for(t in 1:length(v0_t)) {
    rate_mat0[t,] = c(Q(v0_t[t], 0, beta))[c(5, 13, 10, 14, 15)]
    
    int_rate_mat0[t,1] = integrated_Q(v0_t[t], 0, beta[1,])
    int_rate_mat0[t,2] = integrated_Q(v0_t[t], 0, beta[2,])
    int_rate_mat0[t,3] = integrated_Q(v0_t[t], 0, beta[3,])
    int_rate_mat0[t,4] = integrated_Q(v0_t[t], 0, beta[4,])
    int_rate_mat0[t,5] = integrated_Q(v0_t[t], 0, beta[5,])
}

# Plotting the results from the packages (rates)
labels = c("S1->S2 (sex = 1)", "S1->S4 (sex = 1)", "S2->S3 (sex = 1)",
           "S2->S4 (sex = 1)", "S3->S4 (sex = 1)",
           "S1->S2 (sex = 0)", "S1->S4 (sex = 0)", "S2->S3 (sex = 0)",
           "S2->S4 (sex = 0)", "S3->S4 (sex = 0)")

pdf("Plots/rate_plot_AJ.pdf")
par(mfrow=c(3,2))
for(i in 1:5) {
    ylimit1 = c(min(c(v1[[i]], rate_mat1[,i], int_rate_mat1[,i])),
               max(c(v1[[i]], rate_mat1[,i], int_rate_mat1[,i])))
    plot(v1_t, v1[[i]], type = "l", lty = 2, lwd = 2, 
         xlab = "", ylab = "", main = labels[i], col = "red", ylim = ylimit1)
    lines(v1_t, rate_mat1[,i], lty = 5, col = "red")
    lines(v1_t, int_rate_mat1[,i], lty = 1, col = "black")
    if(i == 1) {
        legend(x = "topleft",          # Position
               legend = c("AJ", "rate", "int. rate"),  # Legend texts
               lty = c(2, 5, 1),           # Line types
               col = c("red", "red", "black"),           # Line colors
               lwd = c(2, 1, 1))
    }

    ylimit0 = c(min(c(v0[[i]], rate_mat0[,i], int_rate_mat0[,i])),
                max(c(v0[[i]], rate_mat0[,i], int_rate_mat0[,i])))
    plot(v0_t, v0[[i]], type = "l", lty = 2, lwd = 2, 
         xlab = "", ylab = "", main = labels[i+5], col = "red", ylim = ylimit0)
    lines(v0_t, rate_mat0[,i], lty = 5, col = "red")
    lines(v0_t, int_rate_mat0[,i], lty = 1, col = "black")
    
}
dev.off()

# Plotting the results from the packages (probabilities)
labels_p = c("S1 (sex = 1)", "S2 (sex = 1)", "S3 (sex = 1)",
           "S4 (sex = 1)", 
           "S1 (sex = 0)", "S2 (sex = 0)", "S3 (sex = 0)",
           "S4 (sex = 0)")

pdf("Plots/prob_plot_AJ.pdf")
par(mfrow=c(3,2))
for(i in 1:4) {
    ylimit1 = c(0,1)
    plot(v1_t, p1[[i]], type = "l", lty = 2, lwd = 2, 
         xlab = "", ylab = "", main = labels_p[i], col = "red", ylim = ylimit1)
    lines(seq(0, max(v1_t), by = 0.01), P1[[i]], lty = 1, col = "black")
    if(i == 1) {
        legend(x = "topleft",          # Position
               legend = c("AJ", "true"),  # Legend texts
               lty = c(2,1),           # Line types
               col = c("red", "black"),           # Line colors
               lwd = c(2, 1))
    }
    
    ylimit0 = c(0,1)
    plot(v0_t, p0[[i]], type = "l", lty = 2, lwd = 2, 
         xlab = "", ylab = "", main = labels_p[i+4], col = "red", ylim = ylimit0)
    lines(seq(0, max(v0_t), by = 0.01), P0[[i]], lty = 1, col = "black")
    
}
dev.off()

# ------------------------------------------------------------------------------
# Testing if we can retrieve the rates (not integrated rates) ------------------
# ------------------------------------------------------------------------------
q1 = list()
q1[[1]] <- cbind(diff(v1[[1]]), v1_t[-1])
q1[[2]] <- cbind(diff(v1[[2]]), v1_t[-1])
q1[[3]] <- cbind(diff(v1[[3]]), v1_t[-1])
q1[[4]] <- cbind(diff(v1[[4]]), v1_t[-1])
q1[[5]] <- cbind(diff(v1[[5]]), v1_t[-1])

q0 = list()
q0[[1]] <- cbind(diff(v0[[1]]), v0_t[-1])
q0[[2]] <- cbind(diff(v0[[2]]), v0_t[-1])
q0[[3]] <- cbind(diff(v0[[3]]), v0_t[-1])
q0[[4]] <- cbind(diff(v0[[4]]), v0_t[-1])
q0[[5]] <- cbind(diff(v0[[5]]), v0_t[-1])

for(i in 1:5) {
    q1[[i]] = q1[[i]][q1[[i]][,1] != 0, ]
    q0[[i]] = q0[[i]][q0[[i]][,1] != 0, ]
}

r1 = list()
r0 = list()
rate_mat_ind = c(5, 13, 10, 14, 15)
for(i in 1:5) {
    r1_times = q1[[i]][,2]
    r0_times = q0[[i]][,2]
    
    r1[[i]] = matrix(nrow = nrow(q1[[i]]), ncol = 2); r1[[i]][,2] = r1_times
    r0[[i]] = matrix(nrow = nrow(q0[[i]]), ncol = 2); r0[[i]][,2] = r0_times
    
    for(t in 1:length(r1_times)) {
        r1[[i]][t,1] = c(Q(r1_times[t], 1, beta))[rate_mat_ind[i]]
    }
    for(t in 1:length(r0_times)) {
        r0[[i]][t,1] = c(Q(r0_times[t], 0, beta))[rate_mat_ind[i]]
    }
}

pdf("Plots/rate_plot_AJ_real.pdf")
par(mfrow=c(3,2))
for(i in 1:5) {
    plot(q1[[i]][,2], log(q1[[i]][,1]), type = "l", lwd = 1, 
         xlab = "", ylab = "", main = labels[i], col = "red")
    lines(r1[[i]][,2], log(r1[[i]][,1]), lty = 1, col = "black")
    if(i == 1) {
        legend(x = "topleft",          # Position
               legend = c("AJ", "true"),  # Legend texts
               lty = c(1,1),           # Line types
               col = c("red", "black"),           # Line colors
               lwd = c(1, 1))
    }
    
    plot(q0[[i]][,2], log(q0[[i]][,1]), type = "l", lwd = 1, 
         xlab = "", ylab = "", main = labels[i+5], col = "red")
    lines(r0[[i]][,2], log(r0[[i]][,1]), lty = 1, col = "black")
    
}
dev.off()


v_both = list()
for(jj in 1:5) {
    df0 = data.frame("y" = log(q0[[jj]][,1]), "t" = q0[[jj]][,2], 
                     "x" = rep(0, nrow(q0[[jj]])))
    df1 = data.frame("y" = log(q1[[jj]][,1]), "t" = q1[[jj]][,2], 
                     "x" = rep(1, nrow(q1[[jj]])))
    
    df_tot = rbind(df0, df1)
    
    df_tot = df_tot[is.finite(df_tot$y), ]

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

# ------------------------------------------------------------------------------
# Split data based on covariate value ------------------------------------------
# ------------------------------------------------------------------------------
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

labels = c("S1->S2 (sex = 1)", "S1->S4 (sex = 1)", "S2->S3 (sex = 1)",
           "S2->S4 (sex = 1)", "S3->S4 (sex = 1)",
           "S1->S2 (sex = 0)", "S1->S4 (sex = 0)", "S2->S3 (sex = 0)",
           "S2->S4 (sex = 0)", "S3->S4 (sex = 0)")

pdf("Plots/compare_split_rate.pdf")
par(mfrow=c(3,2))
for(i in 1:5) {
    ylimit1 = c(min(c(v1[[i]], v1_split[[i]])),
                max(c(v1[[i]], v1_split[[i]])))
    plot(v1_t, v1[[i]], type = "l", lwd = 3,
         xlab = "", ylab = "", main = labels[i], col = "red", ylim = ylimit1)
    lines(v1_t_split, v1_split[[i]], lty = 1, col = "black")
    
    
    ylimit0 = c(min(c(v0[[i]], v0_split[[i]])),
                max(c(v0[[i]], v0_split[[i]])))
    plot(v0_t, v0[[i]], type = "l", lwd = 3,
         xlab = "", ylab = "", main = labels[i+5], col = "red", ylim = ylimit0)
    lines(v0_t_split, v0_split[[i]], lty = 1, col = "black")
    
}
dev.off()

# ------------------------------------------------------------------------------
# Using optim to get parameter estimates ---------------------------------------
# ------------------------------------------------------------------------------

min_residuals <- function(data, par) {
    with(data, sum(( (1/par[2]) * (exp(par[1] + par[2]*t + par[3]*x) - exp(par[1] + par[3]*x)) - y)^2))
}

optim_coeff = matrix(0, nrow = 5, ncol = 3)
for(i in 1:5) {
    
    df0 = data.frame("y" = v0[[i]], "t" = v0_t, 
                     "x" = rep(0, length(v0_t)))
    df1 = data.frame("y" = v1[[i]], "t" = v1_t, 
                     "x" = rep(1, length(v1_t)))
    
    df_tot = rbind(df0, df1)
    
    init_par = beta[i,]
    
    op = optim(par=init_par, fn=min_residuals, data=df_tot)
    
    if(op$convergence != 0) {print("no convergence")}
    
    optim_coeff[i, ] = op$par
}

optim_coeff_split = matrix(0, nrow = 5, ncol = 3)
for(i in 1:5) {
    
    df0 = data.frame("y" = v0_split[[i]], "t" = v0_t_split, 
                     "x" = rep(0, length(v0_t_split)))
    df1 = data.frame("y" = v1_split[[i]], "t" = v1_t_split, 
                     "x" = rep(1, length(v1_t_split)))
    
    df_tot = rbind(df0, df1)
    
    init_par = beta[i,]
    
    op = optim(par=init_par, fn=min_residuals, data=df_tot)
    
    if(op$convergence != 0) {print("no convergence")}
    
    optim_coeff_split[i, ] = op$par
}





