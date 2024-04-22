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
init_par = c(c(matrix(c(-2.476429,  0.10511210, -0.4760257,
                        -2.870750, -0.34870395, -0.3987754,
                        -1.020912, -0.15426747, -0.3329449,
                        -3.280449,  0.17704181, -0.2658490,
                        -2.229799,  0.08724047,  0.6220039), ncol=3, byrow=T)),
                  c(-4.782532, -1.108108, -2.621567, -2.269586),
                  c(-6.402050, -8.186261))

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)

# Defining the mean and variance for the flat Gaussian prior
prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

data_files <- c("Year/", "YearTwo/", "Month/")

folder = 1

load(paste0("real_cav_analysis/Data/cavData_year.rda"))

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