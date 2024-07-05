library(latex2exp)
library(tidyverse)
library(gridExtra)

load('real_cav_analysis/Model_out/deSolve/mcmc_out_10.rda')
chain = mcmc_out$chain[10000:25001, ]
ind_keep = seq(1, nrow(chain), by=10)
chain = chain[ind_keep, ]

true_par = colMeans(chain)
# The true values corresponding to the slope coefficient on time are scaled by
# 3 in order to magnify the effect of the different approaches to modelling 
# continuous time HMMs.
true_par[6:10] = 3 * true_par[6:10]

par_index = list(beta=1:15, misclass=16:19, pi_logit=20:21)

labels <- c(TeX(r'($\hat{\beta}_{0,1}:$ Baseline coefficient, transition state 1: "no CAV" $\to$ state 2: "mild CAV")'),
            TeX(r'($\hat{\beta}_{0,2}:$ Baseline coefficient, transition state 1: "no CAV" $\to$ state 4: "death")'),
            TeX(r'($\hat{\beta}_{0,3}:$ Baseline coefficient, transition state 2: "mild CAV" $\to$ state 3: "severe CAV")'),
            TeX(r'($\hat{\beta}_{0,4}:$ Baseline coefficient, transition state 2: "mild CAV" $\to$ state 4: "death")'),
            TeX(r'($\hat{\beta}_{0,5}:$ Baseline coefficient, transition state 3: "severe CAV" $\to$ state 4: "death")'),
            TeX(r'($\hat{\beta}_{1,1}:$ Time coefficient for transition state 1: "no CAV" $\to$ state 2: "mild CAV")'),
            TeX(r'($\hat{\beta}_{1,2}:$ Time coefficient for transition state 1: "no CAV" $\to$ state 4: "death")'),
            TeX(r'($\hat{\beta}_{1,3}:$ Time coefficient for transition state 2: "mild CAV" $\to$ state 3: "severe CAV")'),
            TeX(r'($\hat{\beta}_{1,4}:$ Time coefficient for transition state 2: "mild CAV" $\to$ state 4: "death")'),
            TeX(r'($\hat{\beta}_{1,5}:$ Time coefficient for transition state 3: "severe CAV" $\to$ state 4: "death")'),
            TeX(r'($\hat{\beta}_{2,1}:$ Sex coefficient for transition state 1: "no CAV" $\to$ state 2: "mild CAV")'),
            TeX(r'($\hat{\beta}_{2,2}:$ Sex coefficient for transition state 1: "no CAV" $\to$ state 4: "death")'),
            TeX(r'($\hat{\beta}_{2,3}:$ Sex coefficient for transition state 2: "mild CAV" $\to$ state 3: "severe CAV")'),
            TeX(r'($\hat{\beta}_{2,4}:$ Sex coefficient for transition state 2: "mild CAV" $\to$ state 4: "death")'),
            TeX(r'($\hat{\beta}_{2,5}:$ Sex coefficient for transition state 3: "severe CAV" $\to$ state 4: "death")'),
            TeX(r'(logit $P$(obs. state 2: "mild CAV" $|$ true state 1: "no CAV"))'),
            TeX(r'(logit $P$(obs. state 1: "no CAV" $|$ true state 2: "mild CAV"))'),
            TeX(r'(logit $P$(obs. state 3: "severe CAV" $|$ true state 2: "mild CAV"))'),
            TeX(r'(logit $P$(obs. state 2: "mild CAV" $|$ true state 3: "severe CAV"))'),
            TeX(r'(logit $P$(initial state 2: "mild CAV"))'), TeX(r'(logit $P$(initial state 3: "severe CAV"))'))

# --------------------------------------------------------------
# -------- load and format posterior means and coverage --------
# --------------------------------------------------------------
load("sim_cav_time_inhomog/Plots/cov_df_deSolve.rda")
cov_final = cov_df
load("sim_cav_time_inhomog/Plots/cov_df_nhm.rda")
cov_final = cbind(cov_final, cov_df)
load("sim_cav_time_inhomog/Plots/cov_df_Month.rda")
cov_final = cbind(cov_final, cov_df)
load("sim_cav_time_inhomog/Plots/cov_df_Year.rda")
cov_final = cbind(cov_final, cov_df)
load("sim_cav_time_inhomog/Plots/cov_df_YearTwo.rda")
cov_final = cbind(cov_final, cov_df)
load("sim_cav_time_inhomog/Plots/cov_df_msm.rda")
cov_final = cbind(cov_final, cov_df)
colnames(cov_final) = c("deSolve", "nhm", "expm_Month", "expm_Year", "expm_BiYear",
                        "msm_Month", "msm_year","msm_BiYear")
cov_final = round(cov_final, digits = 3)

post_means_final = vector(mode = 'list', length = 9)
load("sim_cav_time_inhomog/Plots/post_means_deSolve.rda")
post_means_final[[1]] = post_means
load("sim_cav_time_inhomog/Plots/post_means_nhm.rda")
post_means_final[[2]] = post_means
load("sim_cav_time_inhomog/Plots/post_means_Month.rda")
post_means_final[[3]] = post_means
load("sim_cav_time_inhomog/Plots/post_means_Year.rda")
post_means_final[[4]] = post_means
load("sim_cav_time_inhomog/Plots/post_means_YearTwo.rda")
post_means_final[[5]] = post_means
load("sim_cav_time_inhomog/Plots/post_means_msm.rda")
post_means_final[[6]] = post_means[[1]]; post_means_final[[7]] = post_means[[2]]
post_means_final[[8]] = post_means[[3]]
load("sim_cav_time_inhomog/Plots/post_means_aj.rda")
post_means_final[[9]] = post_means

VP <- vector(mode="list", length = length(labels))
for(r in 1:length(labels)) {

    yVar = disc_type = x_label = NULL

    if(r %in% par_index$beta) {
        yVar = c(post_means_final[[1]][,r], post_means_final[[2]][,r], post_means_final[[3]][,r],
                 post_means_final[[4]][,r], post_means_final[[5]][,r], post_means_final[[6]][,r],
                 post_means_final[[7]][,r], post_means_final[[8]][,r], post_means_final[[9]][,r])
        
        disc_type = c(rep(paste0(TeX(r'((A))'), "\n", cov_final[r,1]), nrow(post_means_final[[1]])),
                      rep(paste0(TeX(r'((B))'), "\n", cov_final[r,2]), nrow(post_means_final[[2]])),
                      rep(paste0(TeX(r'((C): d = 1/6)'), "\n", cov_final[r,3]), nrow(post_means_final[[3]])),
                      rep(paste0(TeX(r'((C): d = 1)'), "\n", cov_final[r,4]), nrow(post_means_final[[4]])),
                      rep(paste0(TeX(r'((C): d = 2)'), "\n",cov_final[r,5]), nrow(post_means_final[[5]])),
                      rep(paste0(TeX(r'((D): d = 1/6)'), "\n",cov_final[r,6]), nrow(post_means_final[[6]])),
                      rep(paste0(TeX(r'((D): d = 1)'), "\n",cov_final[r,7]), nrow(post_means_final[[7]])),
                      rep(paste0(TeX(r'((D): d = 2)'), "\n",cov_final[r,8]), nrow(post_means_final[[8]])),
                      rep(paste0(TeX(r'((E))'), "\n", " "), nrow(post_means_final[[9]])))
        
        Method = c(rep("ODE w/ Post. Means", nrow(post_means_final[[1]])),
                   rep("ODE w/ MLE", nrow(post_means_final[[2]])),
                   rep("Expm w/ Post. Means", nrow(post_means_final[[3]]) + nrow(post_means_final[[4]]) + nrow(post_means_final[[5]])),
                   rep("Expm w/ MLE", nrow(post_means_final[[6]]) + nrow(post_means_final[[7]]) + nrow(post_means_final[[8]])),
                   rep("AJ Estimator", nrow(post_means_final[[9]])))
    } else {
        yVar = c(post_means_final[[1]][,r], post_means_final[[2]][,r], post_means_final[[3]][,r],
                 post_means_final[[4]][,r], post_means_final[[5]][,r], post_means_final[[6]][,r],
                 post_means_final[[7]][,r], post_means_final[[8]][,r])
        
        disc_type = c(rep(paste0(TeX(r'((A))'), "\n", cov_final[r,1]), nrow(post_means_final[[1]])),
                      rep(paste0(TeX(r'((B))'), "\n", cov_final[r,2]), nrow(post_means_final[[2]])),
                      rep(paste0(TeX(r'((C): d = 1/6)'), "\n", cov_final[r,3]), nrow(post_means_final[[3]])),
                      rep(paste0(TeX(r'((C): d = 1)'), "\n", cov_final[r,4]), nrow(post_means_final[[4]])),
                      rep(paste0(TeX(r'((C): d = 2)'), "\n",cov_final[r,5]), nrow(post_means_final[[5]])),
                      rep(paste0(TeX(r'((D): d = 1/6)'), "\n",cov_final[r,6]), nrow(post_means_final[[6]])),
                      rep(paste0(TeX(r'((D): d = 1)'), "\n",cov_final[r,7]), nrow(post_means_final[[7]])),
                      rep(paste0(TeX(r'((D): d = 2)'), "\n",cov_final[r,8]), nrow(post_means_final[[8]])))
        
        Method = c(rep("ODE w/ Post. Means", nrow(post_means_final[[1]])),
                   rep("ODE w/ MLE", nrow(post_means_final[[2]])),
                   rep("Expm w/ Post. Means", nrow(post_means_final[[3]]) + nrow(post_means_final[[4]]) + nrow(post_means_final[[5]])),
                   rep("Expm w/ MLE", nrow(post_means_final[[6]]) + nrow(post_means_final[[7]]) + nrow(post_means_final[[8]])))   
    }

    plot_df = data.frame(yVar = yVar, disc_type = disc_type, Method = Method)
    plot_df$disc_type <- factor(plot_df$disc_type, levels = unique(plot_df$disc_type))
    VP[[r]] = ggplot(plot_df, aes(x=disc_type, y = yVar, fill = Method)) +
      geom_violin(trim=FALSE) +
      geom_boxplot(width=0.1) +
      ggtitle(labels[r]) +
      ylab("") + #paste0("Parameter Value: ", round(true_par[r], 3))
      xlab("") +
      geom_hline(yintercept=true_par[r], linetype="dashed", color = "red", size=1.5) +
      theme(text = element_text(size = 35), legend.position = "none", 
            panel.background = element_rect(fill = "white", colour = "grey", 
                                    linetype = 'solid', size = 1),
            panel.grid.major = element_line(linetype = 'solid',
                                colour = "grey") )
}

png("visualizations/Plots/violinPlots_final_b1.png", width = 1600, height = 700)
print(VP[[1]]);
dev.off()
png("visualizations/Plots/violinPlots_final_b2.png", width = 1600, height = 700)
print(VP[[2]]);
dev.off()
png("visualizations/Plots/violinPlots_final_b3.png", width = 1600, height = 700)
print(VP[[3]]);
dev.off()
png("visualizations/Plots/violinPlots_final_b4.png", width = 1600, height = 700)
print(VP[[4]]);
dev.off()
png("visualizations/Plots/violinPlots_final_b5.png", width = 1600, height = 700)
print(VP[[5]]);
dev.off()

for(i in 6:length(labels)) {
  png(paste0("visualizations/Plots/Supplement/violinPlots_final_",i,".png"), width = 1600, height = 700)
  print(VP[[i]]);
  dev.off()
}