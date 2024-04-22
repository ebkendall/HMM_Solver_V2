# The goal of this R file is to compare the computation time of the likelihood
# evaluation using the numerical ODE solution compared to the three different 
# matrix-exponential-based solutions. The parameter values used for all
# computations will be the posterior means of the last 15,000 steps (thinned to
# keep every 10th row) of seed 10 from running the script 
# real_ecog_analysis/mcmc_runfile.r 10

load('real_ecog_analysis/Model_out/mcmc_out_10.rda')
chain = mcmc_out$chain[10000:25001, ]
ind_keep = seq(1, nrow(chain), by=10)
chain = chain[ind_keep, ]

init_par = colMeans(chain)

library(foreach, quietly=T)
library(doParallel, quietly=T)

par_index = list( beta=1:24, misclass = 25:30, pi_logit=31:33, l_delta = 34:37, 
                  l_theta=38:41, l_alpha=42:45, l_beta=46:49)

prior_mean = c(c(matrix(c(-8, 0,
                          -8, 0,
                          -8, 0,
                          -8, 0,
                          -8, 0,
                          -8, 0,
                          -8, 0,
                          -8, 0,
                          -8, 0,
                          -8, 0,
                          -8, 0,
                          -8, 0), ncol=2, byrow = T)),
                    c(0, 0, 0, 0, 0, 0),
                    c(0, 0, 0),
                    c(0, 0, 0, 0, 
                      0, 0, 0, 0, 
                      0, 0, 0, 0,
                      0, 0, 0, 0))

prior_sd = c(c(matrix(c(2, 20,
                        2, 20,
                        2, 20,
                        2, 20,
                        2, 20,
                        2, 20,
                        2, 20,
                        2, 20,
                        2, 20,
                        2, 20,
                        2, 20,
                        2, 20), ncol=2, byrow = T)), 
                  c(20, 20, 20, 20, 20, 20),
                  c(20, 20, 20),
                  c(20, 20, 20, 20, 
                    20, 20, 20, 20, 
                    20, 20, 20, 20,
                    20, 20, 20, 20))

prior_par = data.frame( prior_mean= prior_mean,
                        prior_sd= prior_sd)

n_cores = 1
cl <- makeCluster(n_cores, outfile="")
registerDoParallel(cl)

e_time = rep(NA, 2)

load('real_ecog_analysis/Data_format/mice_format_sub_total_split.rda')
# Excluding these data because the original EEG data for these indices exhibited 
# behavior not consistent with the rest of the data.
mice_format = mice_format[!(mice_format$ptnum %in% 34:59), ]

# Random subset of 50 mice for the analysis (reduce size for computation efficiency)
sub_ind = c(101,  77,   3,  97,  70,  84, 113,  82,  30,  88,   2,  24,   1,  
             61,  27,  74,  99, 104,  32,  73, 110,  15,  22,  90,  98,  75,  
             80,  11,  16,  14,  21,  78, 102,  91, 103,  96,  13,  20,  12,
             29,   5,  67,  68,  85,   8,  10,  60, 106,  92,  87)
 
mice_format = mice_format[(mice_format$ptnum %in% sub_ind), ]

# Numerical ODE solver --------------------------------------------------------
source('real_ecog_analysis/mcmc_routine.r')

temp_data = as.matrix(mice_format); rownames(temp_data) = NULL
id = temp_data[,"ptnum"]
y_1 = temp_data[,"state"]
y_2 = temp_data[,c("delta", "theta", "alpha", "beta"), drop=F]
t = temp_data[,"t1"]

s_time = Sys.time()
l = fn_log_post_continuous(init_par, prior_par, par_index, y_1, y_2, t, id)
e_time[1] = Sys.time() - s_time
# -----------------------------------------------------------------------------

# Matrix Exponential ----------------------------------------------------------
source('real_ecog_analysis/mcmc_routine_expm.r')

# Setting the discrete time to match the observed times
disc_time = mice_format$t1
mice_format = cbind(mice_format, disc_time)

temp_data = as.matrix(mice_format); rownames(temp_data) = NULL
id = temp_data[,"ptnum"]
y_1 = temp_data[,"state"]
y_2 = temp_data[,c("delta", "theta", "alpha", "beta"), drop=F]
t = temp_data[,c("t1", "disc_time")]

s_time = Sys.time()
l = fn_log_post_continuous(init_par, prior_par, par_index, y_1, y_2, t, id)
e_time[2] = Sys.time() - s_time
# -----------------------------------------------------------------------------

names(e_time) = c('continuous', 'discrete')
print(e_time)

stopCluster(cl)