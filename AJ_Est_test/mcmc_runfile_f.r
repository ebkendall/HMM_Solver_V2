# source("supplement_code/mcmc_routine.r")
source("mcmc_routine.r")

# args = commandArgs(TRUE)
# ind = as.numeric(args[1])
# case_num = as.numeric(args[2])
ind = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
case_num = 2

# case_num = 1 --> not exact transition time nor all transitions observed
# case_num = 2 --> exact transition times
# case_num = 3 --> exact transition times and all transitions observed

set.seed(ind)
print(ind)

# We will start the MCMC at the true values 
# load('real_cav_analysis/Model_out/deSolve/mcmc_out_10.rda')
load('mcmc_out_10.rda')
chain = mcmc_out$chain[10000:25001, ]
ind_keep = seq(1, nrow(chain), by=10)
chain = chain[ind_keep, ]
init_par = colMeans(chain)
init_par[6:10] = 3 * init_par[6:10]
init_par[7] = init_par[8]
init_par[8] = init_par[6]

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)

load(paste0("DataOut/cavData_case", case_num, "_it", ind, ".rda"))

temp_data = as.matrix(cavData); rownames(temp_data) = NULL
id = temp_data[,"ptnum"]
y = temp_data[,"state"]
x = temp_data[, c("disc_time", "sex"), drop=F]
t = temp_data[,"years"]
steps = 20000
burnin = 5000
n_cores = 16
disc = F

# Center time
mean_t = mean(t)
t = t - mean_t
init_par[1:5] = init_par[1:5] + init_par[6:10] * mean_t


# Defining the mean and variance for the flat Gaussian prior
prior_par = data.frame( prior_mean=init_par[c(par_index$beta, par_index$pi_logit)],
                        prior_sd=rep( 100, length(init_par[c(par_index$beta, par_index$pi_logit)])))

s_time = Sys.time()

mcmc_out = mcmc_routine(y, x, t, id, init_par, prior_par, par_index,
                        steps, burnin, n_cores, disc, case_num)

mcmc_out$mean_t = mean_t

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = paste0('Model_out/mcmc_out_case', case_num, "_it", ind, '.rda'))
