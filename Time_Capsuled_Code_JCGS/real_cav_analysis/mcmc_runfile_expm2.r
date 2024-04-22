source("real_cav_analysis/mcmc_routine.r")

args = commandArgs(TRUE)

ind = as.numeric(args[1]) 
set.seed(ind)
print(ind)

# These initial parameters where from running the MCMC for 30,000 steps with
# 5,000 burnin a total of 3 times. After each run of the MCMC routine, one of 
# the 10 seeds was selected, and the next MCMC routine had its initial values 
# set to the posterior means of the last 300 steps. These are the posterior 
# means of the last 300 steps using seed #3.
init_par = c(c(matrix(c(-2.467986,  0.10280604, -0.6567030,
                        -3.147562, -0.14336491,  0.3829111,
                        -1.186460, -0.14364439,  0.2181610,
                        -3.246101,  0.16391031, -1.4640613,
                        -1.880993,  0.05508141,  0.9968800), ncol=3, byrow=T)),
                  c(-4.524266, -1.198862, -2.595439, -1.904722),
                  c(-8.163384, -8.960717))

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)

# Defining the mean and variance for the flat Gaussian prior
prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

data_files <- c("Year/", "YearTwo/", "Month/")

folder = 2

load(paste0("real_cav_analysis/Data/cavData_twoYear.rda"))

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
                             'mcmc_out_',ind, '.rda'))