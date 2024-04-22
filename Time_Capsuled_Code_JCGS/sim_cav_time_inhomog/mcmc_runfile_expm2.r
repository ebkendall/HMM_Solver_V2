source("sim_cav_time_inhomog/mcmc_routine.r")

args = commandArgs(TRUE)

ind = as.numeric(args[1]) 
set.seed(ind)
print(ind)

# We will start the MCMC at the true values 
load('real_cav_analysis/Model_out/deSolve/mcmc_out_10.rda')
chain = mcmc_out$chain[10000:25001, ]
ind_keep = seq(1, nrow(chain), by=10)
chain = chain[ind_keep, ]
init_par = colMeans(chain)
init_par[6:10] = 3 * init_par[6:10]

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)

# Defining the mean and variance for the flat Gaussian prior
prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

data_files <- c("Year/", "YearTwo/", "Month/")

folder = 2

load(paste0('sim_cav_time_inhomog/DataOut/', data_files[folder],'cavData',ind,'.rda'))

temp_data = as.matrix(cavData); rownames(temp_data) = NULL
id = temp_data[,"ptnum"]
y = temp_data[,"state"]
x = temp_data[, c("disc_time", "sex"), drop=F]
t = temp_data[,"years"]
steps = 20000
burnin = 5000
n_cores = 16
disc = T

s_time = Sys.time()

mcmc_out = mcmc_routine(y, x, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores, disc)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = paste0('sim_cav_time_inhomog/Model_out/',data_files[folder],
                             'mcmc_out_',ind,'.rda'))
