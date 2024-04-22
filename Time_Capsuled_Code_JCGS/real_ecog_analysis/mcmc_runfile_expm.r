source("real_ecog_analysis/mcmc_routine_expm.r")

args = commandArgs(TRUE)

ind = as.numeric(args[1]) 
set.seed(ind)
print(ind)

# These initial parameters come from running the MCMC for several iterations
# and after each run of the MCMC routine, one of the 10 seeds was selected, 
# and the next MCMC routine had its initial values set to the posterior means of
# the last 300 steps of the chain.

init_par = c(c(matrix(c(4.4190728,   0.1593122,
                        0.6550658,   0.3268463,
                       -9.1233968, -10.2544361,
                        4.5131070,   0.2984503,
                       -6.7026947,   2.7106623,
                       -3.9535204,  -4.2646728,
                        1.6923852,   0.6206300,
                      -13.3731155,   3.5219449,
                      -10.1911832,  -2.5826887,
                      -10.5147401,   2.3642762,
                      -10.5910677,  -8.5523652,
                       -6.0113723,  -4.7344687), ncol=2, byrow = T)),
            c(-0.02277076, -6.85184348, -2.27224070, -4.27184991, -5.91699336, -5.53818311),
            c(0.1731581, -1.6327922, -6.8247678),
            c(  2.787542,    3.518268,    1.823414,    1.189779,
                2.604070,    2.515774,    1.696786,   -5.328094,
                1.739428,    1.562293,    0.662239,   -2.837649,
              0.77940107,  0.72884649,  0.04330907, -1.35632967))

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

# Setting the discrete time to match the observed times
disc_time = mice_format$t1
mice_format = cbind(mice_format, disc_time)

temp_data = as.matrix(mice_format); rownames(temp_data) = NULL
id = temp_data[,"ptnum"]
y_1 = temp_data[,"state"]
y_2 = temp_data[,c("delta", "theta", "alpha", "beta"), drop=F]
t = temp_data[,c("t1", "disc_time")]
steps = 30000
burnin = 5000
n_cores = 20

s_time = Sys.time()

mcmc_out = mcmc_routine(y_1, y_2, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores, ind)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = 
              paste0("real_ecog_analysis/Model_out/mcmc_out_", ind, "_expm.rda"))
