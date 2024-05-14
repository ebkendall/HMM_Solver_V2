source("mcmc_routine.r")

# args = commandArgs(TRUE)
# it = as.numeric(args[1])
# exact_time = as.logical(as.numeric(args[2]))
it = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
exact_time = F

set.seed(it)
print(it)


init_par= c(matrix(c(-2, 0.5, -0.4689827,
                     -2, 0.5, -0.2557522,
                     -3, 0.5,  0.1457067,
                     -3, 0.5,  0.8164156,
                     -3, 0.5, -0.5966361,
                     -3, 0.5,  0.7967794), ncol = 3, byrow = T))

par_index = list( beta=1:18)

# Defining the mean and variance for the flat Gaussian prior
prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

if(exact_time) {
    load(paste0('DataOut/exactTime/cavData', it, '.rda'))
} else {
    load(paste0('DataOut/interTime/cavData', it, '.rda'))
}

temp_data = as.matrix(cavData); rownames(temp_data) = NULL
id = temp_data[,"ptnum"]
y = temp_data[,"state"]
x = temp_data[, c("disc_time", "sex"), drop=F]
t = temp_data[,"years"]

# Center Time
mean_time = mean(t)
t = t - mean_time
init_par[1:6] = init_par[1:6] - mean_time * init_par[7:12]


steps = 10000
burnin = 1000
n_cores = 16


s_time = Sys.time()

mcmc_out = mcmc_routine(y, x, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores, exact_time)

mcmc_out$chain[,1:6] = mcmc_out$chain[,1:6] + mean_time * mcmc_out$chain[,7:12]

e_time = Sys.time() - s_time; print(e_time)

if(exact_time) {
    save(mcmc_out, file = paste0("Model_out/exactTime/mcmc_out_", it, ".rda"))
} else {
    save(mcmc_out, file = paste0("Model_out/interTime/mcmc_out_", it, ".rda"))
}
