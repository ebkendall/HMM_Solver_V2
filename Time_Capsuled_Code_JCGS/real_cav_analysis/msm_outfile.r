library(tidyverse)
library(gridExtra)

args = commandArgs(TRUE)

seed_ind = as.numeric(args[1]) 
print(seed_ind)

month_data = matrix(data=-1, nrow = 1, ncol = 21)
year_data = matrix(data=-1, nrow = 1, ncol = 21)
year_2_data = matrix(data=-1, nrow = 1, ncol = 21)

parEst_msm_MONTH = upper_msm_MONTH = lower_msm_MONTH = 
    parEst_msm_YEAR = upper_msm_YEAR = lower_msm_YEAR = 
    parEst_msm_YEAR2 = upper_msm_YEAR2 = lower_msm_YEAR2 = NULL


i = 1
# ---------------------------------------------------------------------------
load(paste0("real_cav_analysis/Model_out/Month_msm/Output_msm", seed_ind, ".rda"))
month_data[i,] = c(Output_msm$opt$par)

# saving the 95% confidence interval
temp <- NULL
temp_upper_msm <- NULL
temp_lower_msm <- NULL
for(l in 1:3){  
  for(r in c(5,13,10,14,15)){  
    temp <- c( temp, Output_msm$Qmatrices[[l]][r])  
    temp_upper_msm <- c( temp_upper_msm, 
                Output_msm$Qmatrices[[l]][r] + 1.96*Output_msm$QmatricesSE[[l]][r])
    temp_lower_msm <- c( temp_lower_msm, 
                Output_msm$Qmatrices[[l]][r] - 1.96*Output_msm$QmatricesSE[[l]][r])
  }  
}

# Putting misclassification on the logit scale
misclass_vals = Output_msm$Ematrices[[1]][c(5,2,10,7)]
lowerBounds = Output_msm$EmatricesL[[1]][c(5,2,10,7)]
upperBounds = Output_msm$EmatricesU[[1]][c(5,2,10,7)]

lowerBounds[1] = log(lowerBounds[1] / (1-lowerBounds[1]))
lowerBounds[2:3] = log(lowerBounds[2:3] / (1-sum(lowerBounds[2:3])))
lowerBounds[4] = log(lowerBounds[4] / (1-lowerBounds[4]))

upperBounds[1] = log(upperBounds[1] / (1-upperBounds[1]))
upperBounds[2:3] = log(upperBounds[2:3] / (1-sum(upperBounds[2:3])))
upperBounds[4] = log(upperBounds[4] / (1-upperBounds[4]))

temp = c(temp, misclass_vals)
temp_upper_msm <- c( temp_upper_msm, upperBounds)
temp_lower_msm <- c( temp_lower_msm, lowerBounds)

# Putting initial state probabilities on the logit scale
init_vals = Output_msm$opt$par[20:21]
l_bound_init = c(Output_msm$ci[36,1], Output_msm$ci[37,1])
u_bound_init = c(Output_msm$ci[36,2], Output_msm$ci[37,2])

l_bound_init = log(l_bound_init / (1 - sum(l_bound_init)))
u_bound_init = log(u_bound_init / (1 - sum(u_bound_init)))

temp <- c( temp, init_vals)
temp_upper_msm <- c( temp_upper_msm, u_bound_init)
temp_lower_msm <- c( temp_lower_msm, l_bound_init)

parEst_msm_MONTH <- rbind( parEst_msm_MONTH, temp)
upper_msm_MONTH <- rbind( upper_msm_MONTH, temp_upper_msm)
lower_msm_MONTH <- rbind( lower_msm_MONTH, temp_lower_msm)

# ---------------------------------------------------------------------------
load(paste0("real_cav_analysis/Model_out/Year_msm/Output_msm", seed_ind, ".rda"))
year_data[i,] = c(Output_msm$opt$par)

# saving the 95% confidence interval
temp <- NULL
temp_upper_msm <- NULL
temp_lower_msm <- NULL
for(l in 1:3){  
  for(r in c(5,13,10,14,15)){  
    temp <- c( temp, Output_msm$Qmatrices[[l]][r])  
    temp_upper_msm <- c( temp_upper_msm, 
                Output_msm$Qmatrices[[l]][r] + 1.96*Output_msm$QmatricesSE[[l]][r])
    temp_lower_msm <- c( temp_lower_msm, 
                Output_msm$Qmatrices[[l]][r] - 1.96*Output_msm$QmatricesSE[[l]][r])
  }  
}
# Putting misclassification on the logit scale
misclass_vals = Output_msm$Ematrices[[1]][c(5,2,10,7)]
lowerBounds = Output_msm$EmatricesL[[1]][c(5,2,10,7)]
upperBounds = Output_msm$EmatricesU[[1]][c(5,2,10,7)]

lowerBounds[1] = log(lowerBounds[1] / (1-lowerBounds[1]))
lowerBounds[2:3] = log(lowerBounds[2:3] / (1-sum(lowerBounds[2:3])))
lowerBounds[4] = log(lowerBounds[4] / (1-lowerBounds[4]))

upperBounds[1] = log(upperBounds[1] / (1-upperBounds[1]))
upperBounds[2:3] = log(upperBounds[2:3] / (1-sum(upperBounds[2:3])))
upperBounds[4] = log(upperBounds[4] / (1-upperBounds[4]))

temp = c(temp, misclass_vals)
temp_upper_msm <- c( temp_upper_msm, upperBounds)
temp_lower_msm <- c( temp_lower_msm, lowerBounds)

# Putting initial state probabilities on the logit scale
init_vals = Output_msm$opt$par[20:21]
l_bound_init = c(Output_msm$ci[36,1], Output_msm$ci[37,1])
u_bound_init = c(Output_msm$ci[36,2], Output_msm$ci[37,2])

l_bound_init = log(l_bound_init / (1 - sum(l_bound_init)))
u_bound_init = log(u_bound_init / (1 - sum(u_bound_init)))

temp <- c( temp, init_vals)
temp_upper_msm <- c( temp_upper_msm, u_bound_init)
temp_lower_msm <- c( temp_lower_msm, l_bound_init)

parEst_msm_YEAR <- rbind( parEst_msm_YEAR, temp)
upper_msm_YEAR  <- rbind( upper_msm_YEAR, temp_upper_msm)
lower_msm_YEAR  <- rbind( lower_msm_YEAR, temp_lower_msm)

# ---------------------------------------------------------------------------
load(paste0("real_cav_analysis/Model_out/YearTwo_msm/Output_msm", seed_ind, ".rda"))
year_2_data[i,] = c(Output_msm$opt$par)

# saving the 95% confidence interval
temp <- NULL
temp_upper_msm <- NULL
temp_lower_msm <- NULL
for(l in 1:3){  
  for(r in c(5,13,10,14,15)){  
    temp <- c( temp, Output_msm$Qmatrices[[l]][r])  
    temp_upper_msm <- c( temp_upper_msm, 
                Output_msm$Qmatrices[[l]][r] + 1.96*Output_msm$QmatricesSE[[l]][r])
    temp_lower_msm <- c( temp_lower_msm, 
                Output_msm$Qmatrices[[l]][r] - 1.96*Output_msm$QmatricesSE[[l]][r])
  }  
}
# Putting misclassification on the logit scale
misclass_vals = Output_msm$Ematrices[[1]][c(5,2,10,7)]
lowerBounds = Output_msm$EmatricesL[[1]][c(5,2,10,7)]
upperBounds = Output_msm$EmatricesU[[1]][c(5,2,10,7)]

lowerBounds[1] = log(lowerBounds[1] / (1-lowerBounds[1]))
lowerBounds[2:3] = log(lowerBounds[2:3] / (1-sum(lowerBounds[2:3])))
lowerBounds[4] = log(lowerBounds[4] / (1-lowerBounds[4]))

upperBounds[1] = log(upperBounds[1] / (1-upperBounds[1]))
upperBounds[2:3] = log(upperBounds[2:3] / (1-sum(upperBounds[2:3])))
upperBounds[4] = log(upperBounds[4] / (1-upperBounds[4]))

temp = c(temp, misclass_vals)
temp_upper_msm <- c( temp_upper_msm, upperBounds)
temp_lower_msm <- c( temp_lower_msm, lowerBounds)

# Putting initial state probabilities on the logit scale
init_vals = Output_msm$opt$par[20:21]
l_bound_init = c(Output_msm$ci[36,1], Output_msm$ci[37,1])
u_bound_init = c(Output_msm$ci[36,2], Output_msm$ci[37,2])

l_bound_init = log(l_bound_init / (1 - sum(l_bound_init)))
u_bound_init = log(u_bound_init / (1 - sum(u_bound_init)))

temp <- c( temp, init_vals)
temp_upper_msm <- c( temp_upper_msm, u_bound_init)
temp_lower_msm <- c( temp_lower_msm, l_bound_init)

parEst_msm_YEAR2 <- rbind( parEst_msm_YEAR2, temp)
upper_msm_YEAR2 <- rbind( upper_msm_YEAR2, temp_upper_msm)
lower_msm_YEAR2 <- rbind( lower_msm_YEAR2, temp_lower_msm)

msm_info = list(par = parEst_msm_MONTH, upper = upper_msm_MONTH, lower = lower_msm_MONTH)
save(msm_info, file = "real_cav_analysis/Plots/msm_month.rda")
msm_info = list(par = parEst_msm_YEAR, upper = upper_msm_YEAR, lower = lower_msm_YEAR)
save(msm_info, file = "real_cav_analysis/Plots/msm_year.rda")
msm_info = list(par = parEst_msm_YEAR2, upper = upper_msm_YEAR2, lower = lower_msm_YEAR2)
save(msm_info, file = "real_cav_analysis/Plots/msm_yearTwo.rda")

load_names = c("real_cav_analysis/Plots/msm_month.rda", 
               "real_cav_analysis/Plots/msm_year.rda", 
               "real_cav_analysis/Plots/msm_yearTwo.rda")
save_names = c("real_cav_analysis/Plots/cred_set_cumulative_Month_msm.rda", 
               "real_cav_analysis/Plots/cred_set_cumulative_Year_msm.rda", 
               "real_cav_analysis/Plots/cred_set_cumulative_YearTwo_msm.rda")
for(i in 1:3) {
  load(load_names[i])
  upper = c(msm_info[[2]])
  lower = c(msm_info[[3]])

  cred_set_cumulative = cbind(lower, upper)
  save(cred_set_cumulative, file = save_names[i])
}