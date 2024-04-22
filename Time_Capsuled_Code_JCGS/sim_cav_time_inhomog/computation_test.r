# The goal of this R file is to compare the computation time of the likelihood
# evaluation using the numerical ODE solution compared to the three different 
# matrix-exponential-based solutions. Note that each likelihood calculation
# requires loading a separate data frame because the matrix-exponential-based 
# solution requires the addition of censored rows. The parameter values used for all
# computations will be the posterior means of the last 15,000 steps (thinned to
# keep every 10th row) of seed 10 from running the script 
# real_cav_analysis/mcmc_runfile_deSolve.r

library(foreach, quietly=T)
library(doParallel, quietly=T)

source('sim_cav_time_inhomog/mcmc_routine.r')

load('real_cav_analysis/Model_out/deSolve/mcmc_out_10.rda')
chain = mcmc_out$chain[10000:25001, ]
ind_keep = seq(1, nrow(chain), by=10)
chain = chain[ind_keep, ]
init_par = colMeans(chain)
init_par[6:10] = 3 * init_par[6:10]

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)
prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))
n_cores = 3
cl <- makeCluster(n_cores, outfile="")
registerDoParallel(cl)

e_time_full = matrix(NA, 100, 4)
colnames(e_time_full) = c('continuous', 'year', 'bi-year', 'bi-month')
data_files <- c("Year/", "YearTwo/", "Month/")

for (ind in 1:100) {
    print(ind)
    e_time = rep(NA, 4)

    # Numerical ODE solution ---------------------------------------------------
    load(paste0("sim_cav_time_inhomog/DataOut/Continuous/cavData", ind, ".rda"))
    temp_data = as.matrix(cavData); rownames(temp_data) = NULL
    id = temp_data[,"ptnum"]
    y = temp_data[,"state"]
    x = temp_data[, c("disc_time", "sex"), drop=F]
    t = temp_data[,"years"]
    disc = F

    s_time = Sys.time()
    l = fn_log_post(init_par, prior_par, par_index, x, y, t, id, disc)
    e_time[1] = Sys.time() - s_time

    # --------------------------------------------------------------------------

    # Matrix Exponential (year discretization) ---------------------------------
    folder = 1
    load(paste0('sim_cav_time_inhomog/DataOut/', data_files[folder],'cavData',ind,'.rda'))
    temp_data = as.matrix(cavData); rownames(temp_data) = NULL
    id = temp_data[,"ptnum"]
    y = temp_data[,"state"]
    x = temp_data[, c("disc_time", "sex"), drop=F]
    t = temp_data[,"years"]
    disc = T

    s_time = Sys.time()
    l = fn_log_post(init_par, prior_par, par_index, x, y, t, id, disc)
    e_time[2] = Sys.time() - s_time

    # --------------------------------------------------------------------------

    # Matrix Exponential (bi-yearly discretization) ----------------------------
    folder = 2
    load(paste0('sim_cav_time_inhomog/DataOut/', data_files[folder],'cavData',ind,'.rda'))
    temp_data = as.matrix(cavData); rownames(temp_data) = NULL
    id = temp_data[,"ptnum"]
    y = temp_data[,"state"]
    x = temp_data[, c("disc_time", "sex"), drop=F]
    t = temp_data[,"years"]
    disc = T

    s_time = Sys.time()
    l = fn_log_post(init_par, prior_par, par_index, x, y, t, id, disc)
    e_time[3] = Sys.time() - s_time

    # --------------------------------------------------------------------------

    # Matrix Exponential (bi-monthly discretization) ---------------------------
    folder = 3
    load(paste0('sim_cav_time_inhomog/DataOut/', data_files[folder],'cavData',ind,'.rda'))
    temp_data = as.matrix(cavData); rownames(temp_data) = NULL
    id = temp_data[,"ptnum"]
    y = temp_data[,"state"]
    x = temp_data[, c("disc_time", "sex"), drop=F]
    t = temp_data[,"years"]
    disc = T

    s_time = Sys.time()
    l = fn_log_post(init_par, prior_par, par_index, x, y, t, id, disc)
    e_time[4] = Sys.time() - s_time

    # --------------------------------------------------------------------------

    e_time_full[ind,] = e_time
}

print(colMeans(e_time_full))
save(e_time_full, file = "sim_cav_time_inhomog/Model_out/e_time_full.rda")

stopCluster(cl)