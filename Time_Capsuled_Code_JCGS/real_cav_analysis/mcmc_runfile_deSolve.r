source("real_cav_analysis/mcmc_routine.r")

args = commandArgs(TRUE)

ind = as.numeric(args[1]) 
set.seed(ind)
print(ind)

# These initial parameters where from running the MCMC for 30,000 steps with
# 5,000 burnin a total of 3 times. After each run of the MCMC routine, one of 
# the 10 seeds was selected, and the next MCMC routine had its initial values 
# set to the posterior means of the last 300 steps. These are the posterior 
# means of the last 300 steps using seed #1.
init_par = c(c(matrix(c( -2.348358,  0.1054161, -0.4348307,
                         -1.589180, -4.1241792,  0.7553760,
                         -1.445115, -0.0985043,  0.4890515,
                         -1.960385,  0.0327601,  0.1237441,
                         -2.533004,  0.1166776,  1.4994372), ncol=3, byrow=T)),
            c(-4.864549, -1.000305,  -2.511552,  -2.738977),
            c(-6.536839, -5.683834))

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)

# Defining the mean and variance for the flat Gaussian prior
prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

load(paste0("real_cav_analysis/Data/cavData_cont.rda"))

temp_data = as.matrix(cavData); rownames(temp_data) = NULL
id = temp_data[,"ptnum"]
y = temp_data[,"state"]
x = temp_data[, c("disc_time", "sex"), drop=F]
t = temp_data[,"years"]
steps = 30000
burnin = 5000
n_cores = 16
disc = F


s_time = Sys.time()

mcmc_out = mcmc_routine(y, x, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores, disc)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = paste0("real_cav_analysis/Model_out/deSolve/mcmc_out_", ind, ".rda"))
