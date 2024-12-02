# Load necessary packages
library(ggplot2)
library(MASS)      # For negative binomial models
library(hnp)       # For half-normal plots
library(reshape2)  # For data manipulation
library(dplyr)     # For data manipulation

sink("Output.txt")
data_url <- "https://www.ime.usp.br/~giapaula/nitrofen.txt"
data <- read.table(data_url, header = TRUE)
colnames(data) <- c("dose", "tovos")
#data$dose <- factor(data$dose, levels = c(0, 80, 160, 235, 310))
print(head(data))

ggplot(data, aes(x = dose, y = tovos)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  labs(
    title = "Total Number of Hatched Eggs by Nitrofen Dose",
    x = "Nitrofen Dose (mg/L)",
    y = "Total Number of Hatched Eggs"
  ) +
  theme_classic()


ggsave("boxplot_tovos_vs_dose.png")

model1 <- glm(tovos ~ dose, data = data, family = poisson(link = "log"))
model2 <- glm(tovos ~ dose + I(as.numeric(as.character(dose))^2), data = data, family = poisson(link = "log"))
model3 <- glm(tovos ~ dose, data = data, family = poisson(link = "identity"))
model4 <- glm(tovos ~ dose, data = data, family = poisson(link = "sqrt"))
model5 <- glm(tovos ~ dose + I(as.numeric(as.character(dose))^2), data = data, family = poisson(link = "identity"))
model6 <- glm(tovos ~ dose + I(as.numeric(as.character(dose))^2), data = data, family = poisson(link = "sqrt"))

model_list <- list(model1, model2, model3, model4, model5, model6)
model_names <- c("model1", "model2", "model3", "model4", "model5", "model6")
aic_values <- sapply(model_list, AIC)
aic_table <- data.frame(Model = model_names, AIC = aic_values)
print(aic_table)

for (i in 1:length(model_list)) {
  model <- model_list[[i]]
  model_name <- model_names[i]
  
  # Generate hnp plot
  hnp_result <- hnp(
    model,
    how.many.out = TRUE,
    paint.out = TRUE,
    sim = 1000,
    main = paste("Half-Normal Plot of Residuals -", model_name)
  )
  
  # Save the hnp plot
  png_filename <- paste0("hnp_plot_", model_name, ".png")
  png(png_filename)
  plot(hnp_result)
  dev.off()
}

for (i in 1:length(model_list)) {
  model <- model_list[[i]]
  model_name <- model_names[i]
  
  # Compute Cook's Distance
  cooksd <- cooks.distance(model)
  
  # Plot Cook's Distance
  png_filename <- paste0("cooks_distance_", model_name, ".png")
  png(png_filename, width = 800, height = 600)
  plot(cooksd, type = "h", main = paste("Cook's Distance -", model_name), ylab = "Cook's Distance")
  abline(h = 4 / nrow(data), col = "red", lty = 2)
  dev.off()
  
  # Identify Influential Observations
  influential_obs <- which(cooksd > (4 / nrow(data)))
  cat("Model:", model_name, "\n")
  cat("Influential Observations based on Cook's Distance:\n")
  print(influential_obs)
  cat("\n")
}

#Model 2 is chosen to be the best
cat("Before refitting without influential points", "\n")

exp_coefficients <- exp(coef(model2))
exp_confint <- exp(confint(model2))

results_table <- data.frame(
  Coefficient = names(coef(model2)),
  Estimate = coef(model2),
  exp_Estimate = exp_coefficients,
  Lower_CI = confint(model2)[,1],
  Upper_CI = confint(model2)[,2],
  exp_Lower_CI = exp(confint(model2)[,1]),
  exp_Upper_CI = exp(confint(model2)[,2])
)

print(results_table)

data$predicted_counts <- predict(model2, type = "response")

ggplot(data, aes(x = dose, y = tovos)) +
  geom_point(color = "blue", size = 3) +
  geom_line(aes(y = predicted_counts), color = "red", size = 1) +
  labs(
    title = "Observed and Predicted Counts",
    x = "Nitrofen Dose (mg/L)",
    y = "Total Number of Hatched Eggs"
  ) +
  theme_bw()

ggsave("observed_vs_predicted.png")

model <- model_list[[2]]
model_name <- model_names[2]
  
cooksd <- cooks.distance(model)
influential_obs <- which(cooksd > (4 / nrow(data)))

cat("After refitting without influential points", "\n")
updated_data <- data[-influential_obs, ]
# Function to refit model excluding influential observations
refit_model <- function(updated_data) {
  glm.nb(
    formula(model2),
    data = updated_data
  )
}

# Refit the model without influential observations
if (length(influential_obs) > 0) {
  selected_model_no_influential <- refit_model(updated_data)
  
  # Extract coefficients
  original_betas <- coef(model2)
  betas_no_influential <- coef(selected_model_no_influential)
  
  # Find common coefficients
  common_coefficients <- intersect(names(original_betas), names(betas_no_influential))
  
  # Create comparison data frame
  beta_comparison <- data.frame(
    Coefficient = common_coefficients,
    Original = original_betas[common_coefficients],
    No_Influential = betas_no_influential[common_coefficients]
  )
  
  beta_comparison$Abs_Diff <- abs(beta_comparison$Original - beta_comparison$No_Influential)
  beta_comparison$Rel_Abs_Diff <- beta_comparison$Abs_Diff / abs(beta_comparison$Original)
  
  print(beta_comparison)
  
  # Plotting the comparison of coefficients
  beta_long <- melt(beta_comparison, id.vars = "Coefficient",
                    measure.vars = c("Original", "No_Influential"))
  
  png("beta_estimates_comparison.png", width = 800, height = 600)
  p_beta <- ggplot(beta_long, aes(x = Coefficient, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    coord_flip() +
    labs(title = "Comparison of Beta Estimates (Common Coefficients)",
         y = "Beta Value",
         fill = "Model") +
    theme_bw()
  print(p_beta)
  dev.off()
  
  # Absolute differences
  beta_abs_long <- melt(beta_comparison, id.vars = "Coefficient",
                        measure.vars = c("Abs_Diff"))
  
  png("beta_absolute_differences.png", width = 800, height = 600)
  p_abs_diff <- ggplot(beta_abs_long, aes(x = Coefficient, y = value)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = "Absolute Differences in Beta Estimates (Common Coefficients)",
         y = "Absolute Difference") +
    theme_bw()
  print(p_abs_diff)
  dev.off()
  
  # Relative absolute differences
  beta_rel_abs_long <- melt(beta_comparison, id.vars = "Coefficient",
                            measure.vars = c("Rel_Abs_Diff"))
  
  png("beta_relative_absolute_differences.png", width = 800, height = 600)
  p_rel_abs_diff <- ggplot(beta_rel_abs_long, aes(x = Coefficient, y = value)) +
    geom_bar(stat = "identity", fill = "darkorange") +
    coord_flip() +
    labs(title = "Relative Absolute Differences in Beta Estimates (Common Coefficients)",
         y = "Relative Absolute Difference") +
    theme_bw()
  print(p_rel_abs_diff)
  dev.off()
  
  # Compare p-values
  p_values_original <- summary(model2)$coefficients[common_coefficients, "Pr(>|z|)"]
  p_values_no_influential <- summary(selected_model_no_influential)$coefficients[common_coefficients, "Pr(>|z|)"]
  
  p_value_comparison <- data.frame(
    Coefficient = common_coefficients,
    Original = p_values_original,
    No_Influential = p_values_no_influential
  )
  
  p_value_long <- melt(p_value_comparison, id.vars = "Coefficient",
                       variable.name = "Model", value.name = "p_value")
  
  p_value_long$minus_log10_pvalue <- -log10(p_value_long$p_value)
  
  # Handle infinite values (p-values of 0)
  max_log10_pvalue <- max(p_value_long$minus_log10_pvalue[is.finite(p_value_long$minus_log10_pvalue)])
  p_value_long$minus_log10_pvalue[is.infinite(p_value_long$minus_log10_pvalue)] <- max_log10_pvalue + 1
  
  # Calculate thresholds
  alpha <- 0.05
  num_tests <- nrow(p_value_comparison)
  bonferroni_threshold <- alpha / num_tests
  minus_log10_bonferroni <- -log10(bonferroni_threshold)
  minus_log10_standard <- -log10(alpha)
  
  # Create thresholds data frame
  thresholds <- data.frame(
    Threshold = c("Bonferroni Threshold", "Standard Threshold"),
    minus_log10_threshold = c(minus_log10_bonferroni, minus_log10_standard)
  )
  
  # Plot p-values with thresholds
  png("p_values_comparison_plot.png", width = 1000, height = 600)
  p_pvalues <- ggplot(p_value_long, aes(x = Coefficient, y = minus_log10_pvalue, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    geom_hline(data = thresholds, aes(yintercept = minus_log10_threshold, linetype = Threshold, color = Threshold), size = 1) +
    labs(
      x = "Coefficient",
      y = expression(-log[10](p~value)),
      title = "Comparison of p-values Across Models (Common Coefficients)",
      fill = "Model",
      linetype = "Threshold",
      color = "Threshold"
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    ) +
    guides(
      linetype = guide_legend(order = 1),
      color = guide_legend(order = 1),
      fill = guide_legend(order = 2)
    )
  print(p_pvalues)
  dev.off()
  
  print(p_value_comparison)

  exp_coefficients <- exp(coef(selected_model_no_influential))
exp_confint <- exp(confint(selected_model_no_influential))

results_table <- data.frame(
  Coefficient = names(coef(selected_model_no_influential)),
  Estimate = coef(selected_model_no_influential),
  exp_Estimate = exp_coefficients,
  Lower_CI = confint(selected_model_no_influential)[,1],
  Upper_CI = confint(selected_model_no_influential)[,2],
  exp_Lower_CI = exp(confint(selected_model_no_influential)[,1]),
  exp_Upper_CI = exp(confint(selected_model_no_influential)[,2])
)

print(results_table)

updated_data$predicted_counts <- predict(selected_model_no_influential, type = "response")

ggplot(updated_data, aes(x = dose, y = tovos)) +
  geom_point(color = "blue", size = 3) +
  geom_line(aes(y = predicted_counts), color = "red", size = 1) +
  labs(
    title = "Observed and Predicted Counts",
    x = "Nitrofen Dose (mg/L)",
    y = "Total Number of Hatched Eggs"
  ) +
  theme_bw()

ggsave("observed_vs_predicted_no_influential.png")
  
} else {
  cat("No influential observations detected based on Cook's Distance.\n")
}


sink()