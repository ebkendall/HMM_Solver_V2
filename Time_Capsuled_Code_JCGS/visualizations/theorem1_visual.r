library(latex2exp)

# Case 1 ---------------------------------------------------------------------
set.seed(2022)
png("visualizations/Plots/case1.png", width = 500, height = 600)

plot(x = NULL, y = NULL, ylim = c(-1,3), xlim = c(0,4), xlab = TeX(r'($t$)'), ylab = TeX(r'($y$)'),
     main = TeX(r'($t_i = g_d(t_i)$)'), xaxs='i', cex.lab=2, cex.axis=2, cex.main=2)
abline(h = 0, lty = 2, lwd = 0.5); abline(v=0, lty = 2, lwd = 0.5)
t = seq(-1, 5, by = 1)
y = -0.5 + t
points(t, y, cex=2)
abline(lm(y~t))

# adding the step function
steps = stepfun(-1:5, -2:5) 
plot(steps, verticals = F, do.points = F, add = T)

new_t = floor(t)
plotting_df = data.frame("y" = y, "t" = t, "new_t" = new_t)
points(new_t, y, pch = 2, cex=2)
abline(lm(y~new_t), lty = 3, lwd = 3)

dev.off()

# Case 2  ---------------------------------------------------------------------
set.seed(2022)
png("visualizations/Plots/case2.png", width = 500, height = 600)

plot(x = NULL, y = NULL, ylim = c(-1,3), xlim = c(0,4), xlab = TeX(r'($t$)'), ylab = TeX(r'($y$)'),
     main = TeX(r'($t_i = g_d(t_i) + c$)'), xaxs='i', cex.lab=2, cex.axis=2, cex.main=2)
abline(h = 0, lty = 2, lwd = 0.5); abline(v=0, lty = 2, lwd = 0.5)
t = seq(-1.4, 5.6, by = 1)
y = -0.5 + t
points(t, y, cex=2)
abline(lm(y~t))

# adding the step function
steps = stepfun(-1:5, -2:5) 
plot(steps, verticals = F, do.points = F, add = T)

# d-floor for d = 1
new_t = floor(t)
plotting_df = data.frame("y" = y, "t" = t, "new_t" = new_t)
points(new_t, y, pch = 2, cex=2)
abline(lm(y~new_t), lty = 3, lwd = 3)

dev.off()

# Case 3  ---------------------------------------------------------------------
set.seed(2022)
png("visualizations/Plots/case3.png", width = 500, height = 600)

plot(x = NULL, y = NULL, ylim = c(-1,3), xlim = c(0,4), xlab = TeX(r'($t$)'), ylab = TeX(r'($y$)'),
     main = TeX(r'(random $t_i$)'), xaxs='i', cex.lab=2, cex.axis=2, cex.main=2)
abline(h = 0, lty = 2, lwd = 0.5); abline(v=0, lty = 2, lwd = 0.5)

t = sort(runif(100, -2, 6))
y = -0.5 + t
points(t, y, cex=2)
abline(lm(y~t))

# adding the step function
steps = stepfun(-1:5, -2:5) 
plot(steps, verticals = F, do.points = F, add = T)

# d-floor for d = 1
new_t = floor(t)
plotting_df = data.frame("y" = y, "t" = t, "new_t" = new_t)
points(new_t, y, pch = 2, cex=2)
abline(lm(y~new_t), lty = 3, lwd = 3)

dev.off()
