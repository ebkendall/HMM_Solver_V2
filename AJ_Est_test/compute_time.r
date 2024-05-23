# The goal of this R file is to compare the computation time of the likelihood
# evaluation using the numerical ODE solver deSolve vs prodint.
# The parameter values used for all computations will be the posterior means of 
# the last 15,000 steps (thinned to keep every 10th row) of seed 10 from running
# the script mcmc_runfile_deSolve.r

library(foreach, quietly=T)
library(doParallel, quietly=T)
library(AalenJohansen, quietly = T)
library(deSolve, quietly = T)
library(mvtnorm, quietly = T)
library(expm, quietly = T)
library(latex2exp, quietly = T)
library(tidyverse, quietly = T)
library(gridExtra, quietly = T)

# "/Users/ebkendal/Desktop/Research/NCSU/2021/HMM_solver_JASA/Time_Capsuled_code_HMM_data"

# Construct the transition rate matrix
Q <- function(t,x_ik,beta){

    betaMat = matrix(beta, ncol = 3, byrow = F) 
    q1  = exp( c(1,t,x_ik) %*% betaMat[1,] )  # Transition from state 1 to state 2.
    q2  = exp( c(1,t,x_ik) %*% betaMat[2,] )  # Transition from state 1 to death.
    q3  = exp( c(1,t,x_ik) %*% betaMat[3,] )  # Transition from state 2 to state 3.
    q4  = exp( c(1,t,x_ik) %*% betaMat[4,] )  # Transition from state 2 to death.
    q5  = exp( c(1,t,x_ik) %*% betaMat[5,] )  # Transition from state 3 to death.

    qmat = matrix(c( 0,q1, 0,q2,
                    0, 0,q3,q4,
                    0, 0, 0,q5,
                    0, 0, 0, 0),nrow=4,byrow=TRUE)
    diag(qmat) = -rowSums(qmat)

    return(qmat)
}

# Setting up the differential equations for deSolve
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

# Evaluating the log posterior
fn_log_post <- function(pars, prior_par, par_index, x, y, t, id, disc, p_int, step_div) {

    # Initial state probabilities
    init_logit = c( 1, exp(pars[par_index$pi_logit][1]), exp(pars[par_index$pi_logit][2]), 0)
    init = init_logit / sum(init_logit)

    # Misclassification response matrix
    resp_fnc = matrix(c(1, exp(pars[par_index$misclass][1]), 0, 0,
                        exp(pars[par_index$misclass][2]), 1, exp(pars[par_index$misclass][3]), 0,
                        0, exp(pars[par_index$misclass][4]),1, 0,
                        0,   0,   0, 1), ncol=4, byrow=TRUE)
    resp_fnc = resp_fnc / rowSums(resp_fnc)

    beta <- pars[par_index$beta]
    p_ic <- c(p1=1,p2=0,p3=0,p4=0,p5=1,p6=0,p7=0,p8=1,p9=0) # initial condition for deSolve

    # Parallelized computation of the log-likelihood
    log_total_val = foreach(i=unique(id), .combine='+', 
                            .export = c("model_t", "Q"), 
                            .packages = c("deSolve", "expm", "AalenJohansen")) %dopar% {

    val = 1; disc_t_i = NULL
    y_i = y[id == i]                # the observed state
    x_i = x[id == i,"sex",drop = F] # only the sex covariate
    t_i = t[id == i]                # continuous time

    # If time is discretized for matrix exponential (disc = TRUE)
    if(disc) { disc_t_i = x[id == i,"disc_time",drop = F] }

    f_i = init %*% diag(resp_fnc[, y_i[1]])
    log_norm = 0

        for(k in 2:length(t_i)) {

            if(disc) {
                P = expm((t_i[k] - t_i[k-1]) * Q(disc_t_i[k-1], x_i[k-1,], beta))
            } else {
                if(p_int) {
                    t1 = t_i[k-1]
                    t2 = t_i[k]
                    step_size = (t2-t1)/step_div
                    prodInt_list <- prodint(t1, t2, step_size, function(t){Q(t, x_ik = x_i[k,], beta = beta)})
                    
                    P = prodInt_list[[length(prodInt_list)]]
                } else {
                    out <- deSolve::ode(p_ic, times = t_i[(k-1):k], func = model_t,
                                        parms = list(b=beta, x_ik = x_i[k,]))
                    
                    P <- matrix(c(out[2,"p1"], out[2,"p2"], out[2,"p3"], out[2,"p4"],
                                0, out[2,"p5"], out[2,"p6"], out[2,"p7"],
                                0,  0, out[2,"p8"], out[2,"p9"],
                                0,  0,  0,  1), nrow = 4, byrow = T)   
                }
            }

            # Checking if death (state = 4) has occurred
            if(y_i[k] != 4) {
                # Checking if the row is observed or a censored row
                if(y_i[k] != 99) {
                    val = f_i %*% P %*% diag(resp_fnc[, y_i[k]])
                } else { # censor row (only happens when disc == TRUE)
                    val = f_i %*% P %*% diag(rowSums(resp_fnc[, 1:3]))
                }
            } else { # death is observed
                if(disc){
                    val = f_i %*% P %*% Q(disc_t_i[k], x_i[k,], beta) %*% diag(resp_fnc[, y_i[k]])
                } else{
                    val = f_i %*% P %*% Q(t_i[k], x_i[k,], beta) %*% diag(resp_fnc[, y_i[k]])
                }
            }

            norm_val = sqrt(sum(val^2))
            f_i = val / norm_val
            log_norm = log_norm + log(norm_val)
        }

	    return(log(sum(f_i)) + log_norm)
    }

  mean = prior_par$prior_mean
  sd = diag(prior_par$prior_sd)
  log_prior_dens = dmvnorm( x=pars, mean=mean, sigma=sd, log=T)

  return(log_total_val + log_prior_dens)

}

# ------------------------------------------------------------------------------
# Loading parameter values -----------------------------------------------------
# ------------------------------------------------------------------------------

load('real_cav_analysis/Model_out/deSolve/mcmc_out_10.rda')
chain = mcmc_out$chain[10000:25001, ]
ind_keep = seq(1, nrow(chain), by=10)
chain = chain[ind_keep, ]
par_val = colMeans(chain)
par_val[6:10] = 3 * par_val[6:10]

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)

prior_par = data.frame( prior_mean=rep( 0, length(par_val)),
                        prior_sd=rep( 20, length(par_val)))

# ------------------------------------------------------------------------------
# We compare how the prodint function does with various discretizations --------
# ------------------------------------------------------------------------------
step_div = seq(2, 100, by = 2)

# See how the liklihood differs depending on the discretization
load("sim_cav_time_inhomog/DataOut/Continuous/cavData1.rda")
cavData = cavData[cavData$ptnum %in% 1:100, ]

temp_data = as.matrix(cavData); rownames(temp_data) = NULL
id = temp_data[,"ptnum"]
y = temp_data[,"state"]
x = temp_data[, c("disc_time", "sex"), drop=F]
t = temp_data[,"years"]
disc = F

n_cores = 3
cl <- makeCluster(n_cores, outfile="")
registerDoParallel(cl)

s_time = Sys.time()
like_deSolve = fn_log_post(par_val, prior_par, par_index, x, y, t, id, disc, F, 0)
deSolve_time = Sys.time() - s_time

like_prodint = NULL
prodint_time = NULL
p_int = T
for(s in step_div) {
    print(s)
    
    s_time = Sys.time()
    lp = fn_log_post(par_val, prior_par, par_index, x, y, t, id, disc, T, s)
    prodint_time = c(prodint_time, Sys.time() - s_time)
    
    like_prodint = c(like_prodint, lp)
}

plot_df = data.frame(yVar = prodint_time, xVar = step_div)
time_plot = ggplot(plot_df, aes(x=xVar, y = yVar)) +
    geom_line()+
    geom_point() +
    labs(title = TeX(r'(Computation time using prodint)')) +
    ylab(TeX(r'(computation time (sec.) )')) + 
    xlab(TeX(r'($s$, number of subintervals)')) +
    geom_hline(yintercept=deSolve_time, linetype="dashed", color = "red", linewidth=1.5) +
    theme(text = element_text(size = 30), legend.position = "none", 
          panel.background = element_rect(fill = "white", colour = "grey", 
                                          linetype = 'solid', linewidth = 1),
          panel.grid.major = element_line(linetype = 'solid',
                                          colour = "grey") )

png("visualizations/Plots/prodint_time.png", width = 1600, height = 700)
print(time_plot);
dev.off()

diff_like = (like_prodint - like_deSolve)^2

plot_df = data.frame(yVar = diff_like, xVar = step_div)
like_plot = ggplot(plot_df, aes(x=xVar, y = yVar)) +
    geom_line()+
    geom_point() +
    labs(title = TeX(r'(Squared difference in likelihood between using deSolve and prodint)')) +
    ylab(TeX(r'(squared likelihood diff.)')) + 
    xlab(TeX(r'($s$, number of subintervals)')) +
    theme(text = element_text(size = 30), legend.position = "none", 
          panel.background = element_rect(fill = "white", colour = "grey", 
                                          linetype = 'solid', linewidth = 1),
          panel.grid.major = element_line(linetype = 'solid',
                                          colour = "grey"))
png("visualizations/Plots/prodint_like.png", width = 1600, height = 700)
print(like_plot);
dev.off()


stopCluster(cl)
