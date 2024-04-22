library(tidyverse)
library(latex2exp)

# Figure 4 --------------------------------------------------------------------

# Combining all of the credible sets into 1 data frame
method = c(rep("A", 21),rep("B", 21), rep("C(1/6)", 21), rep("D(1/6)", 21),
                        rep("C(1)", 21), rep("D(1)", 21),
                        rep("C(2)", 21), rep("D(2)", 21))
baseline_num = rep(1:21, 8)
plotting_df = data.frame("method" = method, "baseline_num" = baseline_num,
                         "lower" = rep(NA, length(method)), 
                         "upper" = rep(NA, length(method)))
# Loading the credible sets
cred_set = NULL
load('real_cav_analysis/Plots/cred_set_cumulative_deSolve.rda')
cred_set = rbind(cred_set, cred_set_cumulative)
load('real_cav_analysis/Plots/cred_set_cumulative_nhm.rda')
cred_set = rbind(cred_set, cred_set_cumulative)
load('real_cav_analysis/Plots/cred_set_cumulative_Month.rda')
cred_set = rbind(cred_set, cred_set_cumulative)
load('real_cav_analysis/Plots/cred_set_cumulative_Month_msm.rda')
cred_set = rbind(cred_set, cred_set_cumulative)
load('real_cav_analysis/Plots/cred_set_cumulative_Year.rda')
cred_set = rbind(cred_set, cred_set_cumulative)
load('real_cav_analysis/Plots/cred_set_cumulative_Year_msm.rda')
cred_set = rbind(cred_set, cred_set_cumulative)
load('real_cav_analysis/Plots/cred_set_cumulative_YearTwo.rda')
cred_set = rbind(cred_set, cred_set_cumulative)
load('real_cav_analysis/Plots/cred_set_cumulative_YearTwo_msm.rda')
cred_set = rbind(cred_set, cred_set_cumulative)

plotting_df = data.frame("method" = method, "baseline_num" = baseline_num,
                         "lower" = cred_set[,1], "upper" = cred_set[,2])

plotting_df$method <- factor(plotting_df$method, levels = unique(plotting_df$method))
plotting_df$baseline_num <- factor(plotting_df$baseline_num, levels = unique(plotting_df$baseline_num))

# Grouping coefficients
baseline_ind = c(1:5, 22:26, 43:47, 64:68, 85:89, 106:110, 127:131, 148:152)
time_ind = baseline_ind + 5
sex_ind = time_ind + 5
misclass_ind = c(16:19, 37:40, 58:61, 79:82, 100:103, 121:124, 142:145, 163:166)
init_ind = c(20,21, 41,42, 62,63, 83,84, 104,105, 125,126, 146,147, 167,168)

# Baseline --------------------------------------------------------------------
baseline_plot = plotting_df[baseline_ind, ]
rownames(baseline_plot) = NULL
png("visualizations/Plots/baseline_cred_set.png", width = 1600, height = 700)
ggplot(baseline_plot, aes(x = baseline_num, col = method)) + 
    geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.6), size = 1.5, width=0.3) +
    scale_x_discrete(labels = c(TeX(r'($\hat{\beta}_{0,1}:$ 1$\to$2)'), TeX(r'($\hat{\beta}_{0,2}:$ 1$\to$4)'),
                        TeX(r'($\hat{\beta}_{0,3}:$ 2$\to$3)'), TeX(r'($\hat{\beta}_{0,4}:$ 2$\to$4)'),
                        TeX(r'($\hat{\beta}_{0,5}:$ 3$\to$4)'))) +
    labs(title="95% Credible and Confidence Sets", x ="Baseline transition rate parameter", y="") +
    scale_color_discrete(name="Approach", labels = c("(A) ", "(B)", "(C): d = 1/6", "(D): d = 1/6",
                                                 "(C): d = 1", "(D): d = 1",
                                                 "(C): d = 2", "(D): d = 2")) +
    coord_cartesian(clip = "off")+
    theme(text=element_text(size=40),
            axis.text.x = element_text(size=35),
            legend.title=element_text(size=40),
            legend.text=element_text(size=40),
            panel.border = element_blank(),
            plot.background = element_rect( fill = "white"),
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_line(linetype = 'solid',
                                colour = "grey")) +
    geom_vline(xintercept=c(0.5,1.5,2.5,3.5,4.5,5.5), linetype="dotted", size=1.5) 
dev.off()

# time coefficient ------------------------------------------------------------
time_plot = plotting_df[time_ind, ]
rownames(time_plot) = NULL
png("visualizations/Plots/Supplement/time_cred_set.png", width = 1600, height = 700)
ggplot(time_plot, aes(x = baseline_num, col = method)) + 
    geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.6), size = 1.5, width=0.3) +
    scale_x_discrete(labels = c(TeX(r'($\hat{\beta}_{1,1}:$ 1$\to$2)'), TeX(r'($\hat{\beta}_{1,2}:$ 1$\to$4)'),
                        TeX(r'($\hat{\beta}_{1,3}:$ 2$\to$3)'), TeX(r'($\hat{\beta}_{1,4}:$ 2$\to$4)'),
                        TeX(r'($\hat{\beta}_{1,5}:$ 3$\to$4)'))) +
    labs(title="95% Credible and Confidence Sets", x ="Time transition rate parameter", y="") +
    scale_color_discrete(name="Approach", labels = c("(A) ", "(B)", "(C): d = 1/6", "(D): d = 1/6",
                                                 "(C): d = 1", "(D): d = 1",
                                                 "(C): d = 2", "(D): d = 2")) +
    coord_cartesian(clip = "off")+
    theme(text=element_text(size=40),
            axis.text.x = element_text(size=35),
            legend.title=element_text(size=40),
            legend.text=element_text(size=40),
            panel.border = element_blank(),
            plot.background = element_rect( fill = "white"),
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_line(linetype = 'solid',
                                colour = "grey")) +
    geom_vline(xintercept=c(0.5,1.5,2.5,3.5,4.5,5.5), linetype="dotted", size=1.5)
dev.off()

# sex coefficient -------------------------------------------------------------
sex_coeff_plot = plotting_df[sex_ind, ]
rownames(sex_coeff_plot) = NULL
png("visualizations/Plots/Supplement/sex_cred_set.png", width = 1600, height = 700)
ggplot(sex_coeff_plot, aes(x = baseline_num, col = method)) + 
    geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.6), size = 1.5, width=0.3) +
    scale_x_discrete(labels = c(TeX(r'($\hat{\beta}_{2,1}:$ 1$\to$2)'), TeX(r'($\hat{\beta}_{2,2}:$ 1$\to$4)'),
                        TeX(r'($\hat{\beta}_{2,3}:$ 2$\to$3)'), TeX(r'($\hat{\beta}_{2,4}:$ 2$\to$4)'),
                        TeX(r'($\hat{\beta}_{2,5}:$ 3$\to$4)'))) +
    labs(title="95% Credible and Confidence Sets", x ="Sex transition rate parameter", y="") +
    scale_color_discrete(name="Approach", labels = c("(A) ", "(B)", "(C): d = 1/6", "(D): d = 1/6",
                                                 "(C): d = 1", "(D): d = 1",
                                                 "(C): d = 2", "(D): d = 2")) +
    coord_cartesian(clip = "off")+
    theme(text=element_text(size=40),
            axis.text.x = element_text(size=35),
            legend.title=element_text(size=40),
            legend.text=element_text(size=40),
            panel.border = element_blank(),
            plot.background = element_rect( fill = "white"),
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_line(linetype = 'solid',
                                colour = "grey")) +
    geom_vline(xintercept=c(0.5,1.5,2.5,3.5,4.5,5.5), linetype="dotted", size=1.5)
dev.off()

# misclassification -----------------------------------------------------------
misclass_plot = plotting_df[misclass_ind, ]
rownames(misclass_plot) = NULL
png("visualizations/Plots/Supplement/misclass_cred_set.png", width = 1600, height = 700)
ggplot(misclass_plot, aes(x = baseline_num, col = method)) + 
    geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.9), size = 1.5, width=0.3) +
    scale_x_discrete(labels = c(TeX(r'($P$(obs. S2 $|$ true S1))'), TeX(r'($P$(obs. S1 $|$ true S2))'),
                        TeX(r'($P$(obs. S3 $|$ true S2))'), TeX(r'($P$(obs. S2 $|$ true S3))'))) +
    labs(title="95% Credible and Confidence Sets", x ="Logit misclassification probabilities", y="") +
    scale_color_discrete(name="Approach", labels = c("(A) ", "(B)", "(C): d = 1/6", "(D): d = 1/6",
                                                 "(C): d = 1", "(D): d = 1",
                                                 "(C): d = 2", "(D): d = 2")) +
    coord_cartesian(clip = "off")+
    theme(text=element_text(size=40), 
            axis.text.x = element_text(size=25),  
            legend.title=element_text(size=40),
            legend.text=element_text(size=40),
            panel.border = element_blank(),
            plot.background = element_rect( fill = "white"),
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_line(linetype = 'solid',
                                colour = "grey")) +
    geom_vline(xintercept=c(0.5,1.5,2.5,3.5,4.5), linetype="dotted", size=1.5)
dev.off()

# initial state probabilities -------------------------------------------------
init_plot = plotting_df[init_ind, ]
rownames(init_plot) = NULL
png("visualizations/Plots/Supplement/init_cred_set.png", width = 1600, height = 700)
ggplot(init_plot, aes(x = baseline_num, col = method)) + 
    geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.6), size = 1.5, width=0.3) +
    scale_x_discrete(labels = c(TeX(r'($P$(initial S2))'), TeX(r'($P$(initial S3))'))) +
    labs(title="95% Credible and Confidence Sets", x ="Logit initial state probabilities", y="") +
    scale_color_discrete(name="Approach", labels = c("(A) ", "(B)", "(C): d = 1/6", "(D): d = 1/6",
                                                 "(C): d = 1", "(D): d = 1",
                                                 "(C): d = 2", "(D): d = 2")) +
    coord_cartesian(clip = "off")+
    theme(text=element_text(size=40),
            axis.text.x = element_text(size=35),
            legend.title=element_text(size=40),
            legend.text=element_text(size=40),
            panel.border = element_blank(),
            plot.background = element_rect( fill = "white"),
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_line(linetype = 'solid',
                                colour = "grey")) +
    geom_vline(xintercept=c(0.5,1.5,2.5), linetype="dotted", size=1.5)
dev.off()

# Figure 5 --------------------------------------------------------------------

# Combining all of the credible sets into 1 data frame
method = c(rep("(A)", 49), rep("(C): d = .005", 49))
baseline_num = rep(1:49, 2)
plotting_df = data.frame("method" = method, "baseline_num" = baseline_num,
                         "lower" = rep(NA, length(method)), 
                         "upper" = rep(NA, length(method)))
# Loading the credible sets
cred_set = NULL
load('real_ecog_analysis/Plots/cred_set_cumulative_deSolve.rda')
cred_set = rbind(cred_set, cred_set_cumulative)
load('real_ecog_analysis/Plots/cred_set_cumulative_expm.rda')
cred_set = rbind(cred_set, cred_set_cumulative)

plotting_df = data.frame("method" = method, "baseline_num" = baseline_num,
                         "lower" = cred_set[,1], "upper" = cred_set[,2])

plotting_df$method <- factor(plotting_df$method, levels = unique(plotting_df$method))
plotting_df$baseline_num <- factor(plotting_df$baseline_num, levels = unique(plotting_df$baseline_num))

# Grouping coefficients
baseline_ind = c(1,2,4,5,7,8); baseline_ind = c(baseline_ind, baseline_ind + 49)
time_ind = baseline_ind + 12
misclass_ind = 25:30; misclass_ind = c(misclass_ind, misclass_ind + 49)
init_ind  = 31:32; init_ind = c(init_ind, init_ind + 49)
delta_ind = 34:36; delta_ind = c(delta_ind, delta_ind + 49)
theta_ind = 38:40; theta_ind = c(theta_ind, theta_ind + 49)
alpha_ind = 42:44; alpha_ind = c(alpha_ind, alpha_ind + 49)
beta_ind  = 46:48; beta_ind  = c( beta_ind,  beta_ind + 49)

# Baseline --------------------------------------------------------------------
baseline_plot = plotting_df[baseline_ind, ]
rownames(baseline_plot) = NULL
png("visualizations/Plots/baseline_cred_set_ecog.png", width = 1600, height = 700)
ggplot(baseline_plot, aes(x = baseline_num, col = method)) + 
    geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.6), size = 1.5, width=0.3) +
    scale_x_discrete(labels = c(TeX(r'($\hat{\beta}_{0,1}:$ 1$\to$2)'), TeX(r'($\hat{\beta}_{0,2}:$ 1$\to$3)'),
                        TeX(r'($\hat{\beta}_{0,4}:$ 2$\to$1)'), TeX(r'($\hat{\beta}_{0,5}:$ 2$\to$3)'),
                        TeX(r'($\hat{\beta}_{0,7}:$ 3$\to$1)'), TeX(r'($\hat{\beta}_{0,8}:$ 3$\to$2)'))) +
    labs(title="95% Credible Sets", x ="Baseline transition rate parameter", y="") +
    scale_color_discrete(name="Approach", labels = c("(A) ", "(C): d = .005")) +
    coord_cartesian(clip = "off")+
    theme(text=element_text(size=40),
            axis.text.x = element_text(size=35),
            legend.title=element_text(size=40),
            legend.text=element_text(size=40),
            panel.border = element_blank(),
            plot.background = element_rect( fill = "white"),
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_line(linetype = 'solid',
                                colour = "grey")) +
    geom_vline(xintercept=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5), linetype="dotted", size=1.5)
dev.off()

# time --------------------------------------------------------------------
time_plot = plotting_df[time_ind, ]
rownames(time_plot) = NULL
png("visualizations/Plots/time_cred_set_ecog.png", width = 1600, height = 700)
ggplot(time_plot, aes(x = baseline_num, col = method)) + 
    geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.6), size = 1.5, width=0.3) +
    scale_x_discrete(labels = c(TeX(r'($\hat{\beta}_{1,1}:$ 1$\to$2)'), TeX(r'($\hat{\beta}_{1,2}:$ 1$\to$3)'),
                        TeX(r'($\hat{\beta}_{1,4}:$ 2$\to$1)'), TeX(r'($\hat{\beta}_{1,5}:$ 2$\to$3)'),
                        TeX(r'($\hat{\beta}_{1,7}:$ 3$\to$1)'), TeX(r'($\hat{\beta}_{1,8}:$ 3$\to$2)'))) +
    labs(title="95% Credible Sets", x ="Time transition rate parameter", y="") +
    scale_color_discrete(name="Approach", labels = c("(A) ", "(C): d = .005")) +
    coord_cartesian(clip = "off")+
    theme(text=element_text(size=40),
            axis.text.x = element_text(size=35),
            legend.title=element_text(size=40),
            legend.text=element_text(size=40),
            panel.border = element_blank(),
            plot.background = element_rect( fill = "white"),
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_line(linetype = 'solid',
                                colour = "grey")) +
    geom_vline(xintercept=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5), linetype="dotted", size=1.5)
dev.off()

# misclassification -----------------------------------------------------------
misclass_plot = plotting_df[misclass_ind, ]
rownames(misclass_plot) = NULL
png("visualizations/Plots/Supplement/misclass_cred_set_ecog.png", width = 1600, height = 700)
ggplot(misclass_plot, aes(x = baseline_num, col = method)) + 
    geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.8), size = 1.5, width=0.3) +
    scale_x_discrete(labels = c(TeX(r'(obs. S2 $|$ S1)'), TeX(r'(obs. S3 $|$ S1)'),
                                TeX(r'(obs. S1 $|$ S2)'), TeX(r'(obs. S3 $|$ S2)'),
                                TeX(r'(obs. S1 $|$ S3)'), TeX(r'(obs. S2 $|$ S3)'))) +
    labs(title="95% Credible Sets", x ="Logit misclassification probabilities", y="") +
    scale_color_discrete(name="Approach", labels = c("(A) ", "(C): d = .005")) +
    coord_cartesian(clip = "off")+
    theme(text=element_text(size=40),   
            axis.text.x = element_text(size=30),
            legend.title=element_text(size=40),
            legend.text=element_text(size=40),
            panel.border = element_blank(),
            plot.background = element_rect( fill = "white"),
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_line(linetype = 'solid',
                                colour = "grey")) +
    geom_vline(xintercept=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5), linetype="dotted", size=1.5)
dev.off()

# initial state probability ---------------------------------------------------
init_plot = plotting_df[init_ind, ]
rownames(init_plot) = NULL
png("visualizations/Plots/Supplement/init_cred_set_ecog.png", width = 1600, height = 700)
ggplot(init_plot, aes(x = baseline_num, col = method)) + 
    geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.6), size = 1.5, width=0.3) +
    scale_x_discrete(labels = c(TeX(r'($P$(initial state 2))'), TeX(r'($P$(initial state 3))'))) +
    labs(title="95% Credible Sets", x ="Logit initial state probabilities", y="") +
    scale_color_discrete(name="Approach", labels = c("(A) ", "(C): d = .005")) +
    coord_cartesian(clip = "off")+
    theme(text=element_text(size=40),
            axis.text.x = element_text(size=35),
            legend.title=element_text(size=40),
            legend.text=element_text(size=40),
            panel.border = element_blank(),
            plot.background = element_rect( fill = "white"),
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_line(linetype = 'solid',
                                colour = "grey")) +
    geom_vline(xintercept=c(0.5,1.5,2.5), linetype="dotted", size=1.5)
dev.off()

# dirichlet coefficients ---------------------------------------------------
delta_plot = plotting_df[delta_ind, ]
rownames(delta_plot) = NULL
png("visualizations/Plots/Supplement/delta_cred_set_ecog.png", width = 1600, height = 700)
ggplot(delta_plot, aes(x = baseline_num, col = method)) + 
    geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.6), size = 1.5, width=0.3) +
    scale_x_discrete(labels = c(TeX(r'($\hat{\lambda}_{\delta}^{(1)}$ (IS))'), TeX(r'($\hat{\lambda}_{\delta}^{(2)}$ (NREM))'),
                                TeX(r'($\hat{\lambda}_{\delta}^{(3)}$ (REM))'))) +
    labs(title="95% Credible Sets", x =TeX(r'(Dirichlet coefficients for $\delta$ frequency band)'), y="") +
    scale_color_discrete(name="Approach", labels = c("(A) ", "(C): d = .005")) +
    coord_cartesian(clip = "off")+
    theme(text=element_text(size=40),
            axis.text.x = element_text(size=35),
            legend.title=element_text(size=40),
            legend.text=element_text(size=40),
            panel.border = element_blank(),
            plot.background = element_rect( fill = "white"),
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_line(linetype = 'solid',
                                colour = "grey")) +
    geom_vline(xintercept=c(0.5,1.5,2.5,3.5), linetype="dotted", size=1.5)
dev.off()

theta_plot = plotting_df[theta_ind, ]
rownames(theta_plot) = NULL
png("visualizations/Plots/Supplement/theta_cred_set_ecog.png", width = 1600, height = 700)
ggplot(theta_plot, aes(x = baseline_num, col = method)) + 
    geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.6), size = 1.5, width=0.3) +
    scale_x_discrete(labels = c(TeX(r'($\hat{\lambda}_{\theta}^{(1)}$ (IS))'), TeX(r'($\hat{\lambda}_{\theta}^{(2)}$ (NREM))'),
                                TeX(r'($\hat{\lambda}_{\theta}^{(3)}$ (REM))'))) +
    labs(title="95% Credible Sets", x =TeX(r'(Dirichlet coefficients for $\theta$ frequency band)'), y="") +
    scale_color_discrete(name="Approach", labels = c("(A) ", "(C): d = .005")) +
    coord_cartesian(clip = "off")+
    theme(text=element_text(size=40),
            axis.text.x = element_text(size=35),
            legend.title=element_text(size=40),
            legend.text=element_text(size=40),
            panel.border = element_blank(),
            plot.background = element_rect( fill = "white"),
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_line(linetype = 'solid',
                                colour = "grey")) +
    geom_vline(xintercept=c(0.5,1.5,2.5,3.5), linetype="dotted", size=1.5)
dev.off()

alpha_plot = plotting_df[alpha_ind, ]
rownames(alpha_plot) = NULL
png("visualizations/Plots/Supplement/alpha_cred_set_ecog.png", width = 1600, height = 700)
ggplot(alpha_plot, aes(x = baseline_num, col = method)) + 
    geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.6), size = 1.5, width=0.3) +
    scale_x_discrete(labels = c(TeX(r'($\hat{\lambda}_{\alpha}^{(1)}$ (IS))'), TeX(r'($\hat{\lambda}_{\alpha}^{(2)}$ (NREM))'),
                                TeX(r'($\hat{\lambda}_{\alpha}^{(3)}$ (REM))'))) +
    labs(title="95% Credible Sets", x =TeX(r'(Dirichlet coefficients for $\alpha$ frequency band)'), y="") +
    scale_color_discrete(name="Approach", labels = c("(A) ", "(C): d = .005")) +
    coord_cartesian(clip = "off")+
    theme(text=element_text(size=40),
            axis.text.x = element_text(size=35),
            legend.title=element_text(size=40),
            legend.text=element_text(size=40),
            panel.border = element_blank(),
            plot.background = element_rect( fill = "white"),
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_line(linetype = 'solid',
                                colour = "grey")) +
    geom_vline(xintercept=c(0.5,1.5,2.5,3.5), linetype="dotted", size=1.5)
dev.off()

beta_plot = plotting_df[beta_ind, ]
rownames(beta_plot) = NULL
png("visualizations/Plots/Supplement/beta_cred_set_ecog.png", width = 1600, height = 700)
ggplot(beta_plot, aes(x = baseline_num, col = method)) + 
    geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.6), size = 1.5, width=0.3) +
    scale_x_discrete(labels = c(TeX(r'($\hat{\lambda}_{\beta}^{(1)}$ (IS))'), TeX(r'($\hat{\lambda}_{\beta}^{(2)}$ (NREM))'),
                                TeX(r'($\hat{\lambda}_{\beta}^{(3)}$ (REM))'))) +
    labs(title="95% Credible Sets", x =TeX(r'(Dirichlet coefficients for $\beta$ frequency band)'), y="") +
    scale_color_discrete(name="Approach", labels = c("(A) ", "(C): d = .005")) +
    coord_cartesian(clip = "off")+
    theme(text=element_text(size=40),
            axis.text.x = element_text(size=35),
            legend.title=element_text(size=40),
            legend.text=element_text(size=40),
            panel.border = element_blank(),
            plot.background = element_rect( fill = "white"),
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_line(linetype = 'solid',
                                colour = "grey")) +
    geom_vline(xintercept=c(0.5,1.5,2.5,3.5), linetype="dotted", size=1.5)
dev.off()