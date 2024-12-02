library(ggplot2)
library(dplyr)
library(reshape2)

raw_data <- "
2300    65  1   4400   56  0
750   156  1   3000   65  0
4300   100  1   4000   17  0
2600   134  1   1500    7  0
6000    16  1   9000   16  0
10500   108  1   5300   22  0
10000   121  1  10000    3  0
17000     4  1  19000    4  0 
5400    39  1  27000    2  0
7000   143  1  28000    3  0
9400    56  1  31000    8  0
32000    26  1  26000    4  0
35000    22  1  21000    3  0
100000     1  1  79000   30  0
100000     1  1 100000    4  0 
52000     5  1 100000   43  0
100000    65  1
"
sink("gamma_model_output.txt")

data_lines <- unlist(strsplit(raw_data, "\n"))

data_lines <- data_lines[data_lines != ""]

WBC <- numeric()
Tempo <- numeric()
AG <- integer()

for (line in data_lines) {
 
  elements <- unlist(strsplit(line, "\\s+"))
  elements <- elements[elements != ""]
  
  if (length(elements) == 6) {
    WBC <- c(WBC, as.numeric(elements[1]))
    Tempo <- c(Tempo, as.numeric(elements[2]))
    AG <- c(AG, as.integer(elements[3]))
    
    WBC <- c(WBC, as.numeric(elements[4]))
    Tempo <- c(Tempo, as.numeric(elements[5]))
    AG <- c(AG, as.integer(elements[6]))
  } else if (length(elements) == 3) {
    WBC <- c(WBC, as.numeric(elements[1]))
    Tempo <- c(Tempo, as.numeric(elements[2]))
    AG <- c(AG, as.integer(elements[3]))
  } else {
    warning("Unexpected number of elements in line: ", line)
  }
}

data_combined <- data.frame(WBC = WBC, Tempo = Tempo, AG = AG)

data_combined$AG <- factor(data_combined$AG, levels = c(0, 1),
                           labels = c("Negative", "Positive"))

summary(data_combined)

png("histograms_numerical.png")
par(mfrow = c(1, 2))
hist(data_combined$WBC, main = "Histogram of WBC", xlab = "WBC")
hist(data_combined$Tempo, main = "Histogram of Survival Tempo", xlab = "Tempo (weeks)")
dev.off()

png("histogram_log_WBC.png")
hist(log(data_combined$WBC), main = "Histogram of log(WBC)", xlab = "log(WBC)")
dev.off()

boxplot_Tempo_ag <- ggplot(data_combined, aes(x = AG, y = Tempo)) +
  geom_boxplot() +
  labs(title = "Survival Tempo by AG Status", x = "AG Status", y = "Survival Tempo (weeks)") +
  theme_classic()

ggsave("boxplot_Tempo_by_AG.png", plot = boxplot_Tempo_ag, width = 8, height = 6, dpi = 300)

scatter_Tempo_wbc <- ggplot(data_combined, aes(x = log(WBC), y = Tempo, color = AG)) +
  geom_point(size = 3) +
  labs(title = "Survival Tempo vs log(WBC) by AG Status", x = "log(WBC)", y = "Survival Tempo (weeks)") +
  theme_classic()

ggsave("scatter_Tempo_vs_logWBC.png", plot = scatter_Tempo_wbc, width = 8, height = 6, dpi = 300)

gamma_model <- glm(Tempo ~ log(WBC) + AG, family = Gamma(link = "log"), data = data_combined)

summary(gamma_model)


exp_coefficients <- exp(coef(gamma_model))
exp_coefficients

png("gamma_model_diagnostics.png")
par(mfrow = c(2, 2))
plot(gamma_model)
dev.off()

residuals_fitted_plot <- ggplot(data_combined, aes(x = fitted(gamma_model), y = residuals(gamma_model, type = "deviance"))) +
  geom_point() +
  labs(title = "Residuals vs Fitted Values", x = "Fitted Values", y = "Deviance Residuals") +
  theme_classic()

ggsave("residuals_vs_fitted.png", plot = residuals_fitted_plot, width = 8, height = 6, dpi = 300)

png("qqplot_residuals.png")
qqnorm(residuals(gamma_model, type = "deviance"), main = "Normal Q-Q Plot of Deviance Residuals")
qqline(residuals(gamma_model, type = "deviance"))
dev.off()

gamma_dispersion <- sum(residuals(gamma_model, type = "pearson")^2) / gamma_model$df.residual
gamma_dispersion

cat("Gamma Model Dispersion Statistic:", gamma_dispersion, "\n")

cooksd <- cooks.distance(gamma_model)

png("cooks_distance.png")
plot(cooksd, type = "h", main = "Cook's Distance", ylab = "Cook's Distance")
abline(h = 4 / nrow(data_combined), col = "red", lty = 2)
dev.off()

influential_points <- which(cooksd > (4 / nrow(data_combined)))
cat("Influential Points: ", influential_points)

cat("Summary of Gamma GLM:\n")
print(summary(gamma_model))

cat("\nExponentiated Coefficients (Multiplicative Effects):\n")
print(exp_coefficients)

cat("\nInterpretation of Estimates:\n")
cat("The exponentiated coefficients represent the multiplicative change in the mean survival Tempo for a one-unit increase in the predictor.\n")
cat("- For log(WBC): Each one-unit increase in log(WBC) multiplies the mean survival Tempo by", round(exp_coefficients["log(WBC)"], 3), ".\n")
cat("- For AGPositive: Being AG positive multiplies the mean survival Tempo by", round(exp_coefficients["AGPositive"], 3), "compared to AG negative patients.\n")

cat("\nDispersion Statistic:\n")
cat("Gamma Model Dispersion Statistic:", gamma_dispersion, "\n")
if (gamma_dispersion > 1.5) {
  cat("The dispersion statistic is significantly greater than 1, indicating potential overdispersion.\n")
} else {
  cat("The dispersion statistic is close to 1, indicating no significant overdispersion.\n")
}

cat("\nInfluential Observations:\n")
if (length(influential_points) > 0) {
  cat("Observations with high Cook's Distance:\n")
  print(influential_points)
} else {
  cat("No influential observations detected based on Cook's Distance.\n")
}

sink()
sink("influential_observatios_analysis_output.txt")

original_betas <- coef(gamma_model)

data_no_28_32 <- data_combined[-c(28, 32), ]
gamma_model_no_28_32 <- glm(Tempo ~ log(WBC) + AG, family = Gamma(link = "log"), data = data_no_28_32)
betas_no_28_32 <- coef(gamma_model_no_28_32)

data_no_28 <- data_combined[-28, ]
gamma_model_no_28 <- glm(Tempo ~ log(WBC) + AG, family = Gamma(link = "log"), data = data_no_28)
betas_no_28 <- coef(gamma_model_no_28)

data_no_32 <- data_combined[-32, ]
gamma_model_no_32 <- glm(Tempo ~ log(WBC) + AG, family = Gamma(link = "log"), data = data_no_32)
betas_no_32 <- coef(gamma_model_no_32)

beta_comparison <- data.frame(
  Coefficient = names(original_betas),
  Original = original_betas,
  No_28_32 = betas_no_28_32,
  No_28 = betas_no_28,
  No_32 = betas_no_32
)

beta_comparison$Abs_Diff_No_28_32 <- abs(beta_comparison$Original - beta_comparison$No_28_32)
beta_comparison$Abs_Diff_No_28 <- abs(beta_comparison$Original - beta_comparison$No_28)
beta_comparison$Abs_Diff_No_32 <- abs(beta_comparison$Original - beta_comparison$No_32)

beta_comparison$Rel_Abs_Diff_No_28_32 <- beta_comparison$Abs_Diff_No_28_32 / abs(beta_comparison$Original)
beta_comparison$Rel_Abs_Diff_No_28 <- beta_comparison$Abs_Diff_No_28 / abs(beta_comparison$Original)
beta_comparison$Rel_Abs_Diff_No_32 <- beta_comparison$Abs_Diff_No_32 / abs(beta_comparison$Original)

print(beta_comparison)

beta_long <- melt(beta_comparison, id.vars = "Coefficient",
                  measure.vars = c("Original", "No_28_32", "No_28", "No_32"))

p1 <- ggplot(beta_long, aes(x = Coefficient, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  coord_flip() +
  labs(title = "Comparison of Beta Estimates",
       y = "Beta Value",
       fill = "Model") +
  theme_bw()

ggsave(filename = "beta_estimates_comparison.png", plot = p1, width = 8, height = 6, dpi = 300)

beta_abs_long <- melt(beta_comparison, id.vars = "Coefficient",
                      measure.vars = c("Abs_Diff_No_28_32", "Abs_Diff_No_28", "Abs_Diff_No_32"))

p2 <- ggplot(beta_abs_long, aes(x = Coefficient, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  coord_flip() +
  labs(title = "Absolute Differences in Beta Estimates",
       y = "Absolute Difference",
       fill = "Exclusion") +
  theme_bw()

ggsave(filename = "beta_absolute_differences.png", plot = p2, width = 8, height = 6, dpi = 300)

beta_rel_abs_long <- melt(beta_comparison, id.vars = "Coefficient",
                          measure.vars = c("Rel_Abs_Diff_No_28_32", "Rel_Abs_Diff_No_28", "Rel_Abs_Diff_No_32"))

p3 <- ggplot(beta_rel_abs_long, aes(x = Coefficient, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  coord_flip() +
  labs(title = "Relative Absolute Differences in Beta Estimates",
       y = "Relative Absolute Difference",
       fill = "Exclusion") +
  theme_bw()

ggsave(filename = "beta_relative_absolute_differences.png", plot = p3, width = 8, height = 6, dpi = 300)


