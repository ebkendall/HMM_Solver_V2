library(tidyverse)
library(gridExtra)

args = commandArgs(TRUE)

seed_ind = as.numeric(args[1]) 
print(seed_ind)

# ---------------------------------------------------------------------------
load(paste0("real_cav_analysis/Model_out/nhm/nhm_out_", seed_ind, ".rda"))
parEst_nhm = c(gomp_fit$par)

# Undoing the centering by 5 from nhm_runfile.r
parEst_nhm[1:5] = gomp_fit$par[1:5] - 5*gomp_fit$par[6:10]

upper_nhm = lower_nhm = rep(NA, length(parEst_nhm))

# saving the 95% confidence interval
inv_hess = solve(gomp_fit$hess)
for(j in 1:length(parEst_nhm)) {
    upper_nhm[j] = parEst_nhm[j] + 1.96 * sqrt(inv_hess[j,j])
    lower_nhm[j] = parEst_nhm[j] - 1.96 * sqrt(inv_hess[j,j])
}

nhm_info = list(par = parEst_nhm, upper = upper_nhm, lower = lower_nhm)
save(nhm_info, file = "real_cav_analysis/Plots/nhm_info.rda")

cred_set_cumulative = cbind(lower_nhm, upper_nhm)
save(cred_set_cumulative, file = "real_cav_analysis/Plots/cred_set_cumulative_nhm.rda")