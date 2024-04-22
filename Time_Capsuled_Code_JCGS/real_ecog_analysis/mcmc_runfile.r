source("real_ecog_analysis/mcmc_routine.r")

args = commandArgs(TRUE)

ind = as.numeric(args[1]) 
set.seed(ind)
print(ind)

# These initial parameters come from running the MCMC for several iterations
# and after each run of the MCMC routine, one of the 10 seeds was selected, 
# and the next MCMC routine had its initial values set to the posterior means of
# the last 300 steps of the chain.

init_par = c(c(matrix(c(4.3716484,   0.1784532,
                        0.5673361,   0.4719570,
                       -9.2525152, -10.5877664,
                        4.5345324,   0.2765103,
                       -6.6048098,   3.3816701,
                       -4.3269746,  -3.2771459,
                        1.5189832,   0.8938321,
                      -13.1526900,   2.9717948,
                      -11.1898954,  -2.5228124,
                      -11.5263409,   1.7041393,
                      -11.0349643,  -8.9381795,
                       -5.8178577,  -4.1537110), ncol=2, byrow = T)),
            c(-0.01312831, -7.01563467, -2.37405988, -4.33414689, -5.85530074, -5.76750430),
            c(-0.04422176, -1.35732049, -5.86368830),
            c(  2.796010,    3.527022,    1.827246,    1.146568,
                2.604194,    2.517916,    1.707840,   -5.362062,
               1.7372470,   1.5618687,   0.6759148,  -2.7595611,
              0.77855185,  0.73110479,  0.03583184, -1.31704021))

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

temp_data = as.matrix(mice_format); rownames(temp_data) = NULL
id = temp_data[,"ptnum"]
y_1 = temp_data[,"state"]
y_2 = temp_data[,c("delta", "theta", "alpha", "beta"), drop=F]
t = temp_data[,"t1"]
steps = 30000
burnin = 5000
n_cores = 20

s_time = Sys.time()

mcmc_out = mcmc_routine(y_1, y_2, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores, ind)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = paste0("real_ecog_analysis/Model_out/mcmc_out_", ind, ".rda"))
