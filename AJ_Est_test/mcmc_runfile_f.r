source("mcmc_routine.r")

# args = commandArgs(TRUE)
# ind = as.numeric(args[1])
# exact_time = as.logical(as.numeric(args[2]))
ind = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
exact_time = F

set.seed(ind)
print(ind)

# We will start the MCMC at the true values 
load('mcmc_out_10.rda')
chain = mcmc_out$chain[10000:25001, ]
ind_keep = seq(1, nrow(chain), by=10)
chain = chain[ind_keep, ]
init_par = colMeans(chain)
init_par[6:10] = 3 * init_par[6:10]

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)

if(exact_time) {
    load(paste0('DataOut/exactTime/cavData', ind, '.rda'))
} else {
    load(paste0('DataOut/interTime/cavData', ind, '.rda'))
}

temp_data = as.matrix(cavData); rownames(temp_data) = NULL
id = temp_data[,"ptnum"]
y = temp_data[,"state"]
x = temp_data[, c("disc_time", "sex"), drop=F]
t = temp_data[,"years"]
steps = 20000
burnin = 5000
n_cores = 16
disc = F

mean_t = mean(t)
t = t - mean_t
init_par[1:5] = init_par[1:5] + init_par[6:10] * mean_t


# Defining the mean and variance for the flat Gaussian prior
prior_par = data.frame( prior_mean=init_par[c(par_index$beta, par_index$pi_logit)],
                        prior_sd=rep( 100, length(init_par[c(par_index$beta, par_index$pi_logit)])))
# prior_par = data.frame( prior_mean=rep( 0, length(init_par[c(par_index$beta, par_index$pi_logit)])),
#                         prior_sd=rep( 20, length(init_par[c(par_index$beta, par_index$pi_logit)])))

s_time = Sys.time()

mcmc_out = mcmc_routine(y, x, t, id, init_par, prior_par, par_index,
                        steps, burnin, n_cores, disc, exact_time)

mcmc_out$mean_t = mean_t

e_time = Sys.time() - s_time; print(e_time)

if(exact_time) {
    save(mcmc_out, file = paste0("Model_out/exactTime/mcmc_out_", ind, ".rda"))
} else {
    save(mcmc_out, file = paste0("Model_out/interTime/mcmc_out_", ind, ".rda"))
}
