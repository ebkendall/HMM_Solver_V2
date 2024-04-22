source("real_cav_analysis/mcmc_routine.r")

args = commandArgs(TRUE)

ind = as.numeric(args[1]) 
set.seed(ind)
print(ind)

# These initial parameters where from running the MCMC for 30,000 steps with
# 5,000 burnin a total of 3 times. After each run of the MCMC routine, one of 
# the 10 seeds was selected, and the next MCMC routine had its initial values 
# set to the posterior means of the last 300 steps. These are the posterior 
# means of the last 300 steps using seed #6.
init_par = c(c(matrix(c(-2.267967,  0.08672876, -0.48623684,
                        -1.548016, -4.60675016,  0.23368811,
                        -1.518598, -0.10658030,  0.14214568,
                        -2.352224,  0.09542065,  0.04435866,
                        -2.567921,  0.12273417,  1.20522746), ncol=3, byrow=T)),
                  c(-5.3545320, -0.7206255, -2.3635738, -2.0441165),
                  c(-6.441788, -8.570548))

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)

# Defining the mean and variance for the flat Gaussian prior
prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

data_files <- c("Year/", "YearTwo/", "Month/")

folder = 3

load(paste0("real_cav_analysis/Data/cavData_month.rda"))

temp_data = as.matrix(cavData); rownames(temp_data) = NULL
id = temp_data[,"ptnum"]
y = temp_data[,"state"]
x = temp_data[, c("disc_time", "sex"), drop=F]
t = temp_data[,"years"]
steps = 30000
burnin = 5000
n_cores = 16
disc = T

s_time = Sys.time()

mcmc_out = mcmc_routine(y, x, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores, disc)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = paste0('real_cav_analysis/Model_out/',data_files[folder],
                             'mcmc_out_',ind,'.rda'))