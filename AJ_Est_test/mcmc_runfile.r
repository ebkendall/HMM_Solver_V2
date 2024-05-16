source("mcmc_routine.r")

# args = commandArgs(TRUE)
# it = as.numeric(args[1])
# exact_time = as.logical(as.numeric(args[2]))
it = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
exact_time = T

set.seed(it)
print(it)


init_par= c(-2.31617310,  -1.28756312,  -1.10116400,  -2.52367543,  -2.10384797,
              0.27050001, -11.65470594,  -0.49306415,   0.28862090,   0.22731278,
              -0.39079609,  -0.05894252,  -0.32509646,   0.48631653,   0.99565810,
              -5.28923943,  -0.90870027,  -2.40751854,  -2.44696544,  -6.52252202,
              -6.24090500)

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)

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
steps = 10000
burnin = 1000
n_cores = 16


s_time = Sys.time()

mcmc_out = mcmc_routine(y, x, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores)

e_time = Sys.time() - s_time; print(e_time)

if(exact_time) {
    save(mcmc_out, file = paste0("Model_out/exactTime/mcmc_out_", it, ".rda"))
} else {
    save(mcmc_out, file = paste0("Model_out/interTime/mcmc_out_", it, ".rda"))
}
