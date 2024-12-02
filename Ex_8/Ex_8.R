library(MASS)        
library(hnp)         
library(ggplot2)     
library(reshape2) 

sink("output.txt")

data_url <- "https://www.ime.usp.br/~giapaula/canc1.txt"

data_vector <- scan(data_url, what = numeric())

# The data consists of 5 columns:
# 1. Age at first employment (1 to 4)
# 2. Year of first employment (1 to 4)
# 3. Time since first employment (1 to 5)
# 4. Number of cancer cases
# 5. Person-years


data_matrix <- matrix(data_vector, ncol = 5, byrow = TRUE)

data <- as.data.frame(data_matrix)

colnames(data) <- c("AgeFirstEmployment", "YearFirstEmployment",
                    "TimeSinceFirstEmployment", "Cases", "PersonYears")


print(head(data))

# Guessed values, just to make it more realistic
data$AgeFirstEmployment <- factor(data$AgeFirstEmployment,
                                  levels = 1:4,
                                  labels = c("Age<20", "Age20-24", "Age25-29", "Age30+"))

data$YearFirstEmployment <- factor(data$YearFirstEmployment,
                                   levels = 1:4,
                                   labels = c("Before1925", "1925-1929", "1930-1934", "1935+"))

data$TimeSinceFirstEmployment <- factor(data$TimeSinceFirstEmployment,
                                        levels = 1:5,
                                        labels = c("0-9yrs", "10-19yrs", "20-29yrs", "30-39yrs", "40+yrs"))

# Ensure the response variable and offset are numeric
data$Cases <- as.numeric(data$Cases)
data$PersonYears <- as.numeric(data$PersonYears)

# View the structure of the data
print(str(data))

# Fit the initial Poisson regression model with main effects only
main_effects_model <- glm(Cases ~ AgeFirstEmployment + YearFirstEmployment + TimeSinceFirstEmployment,
                          family = poisson(link = "log"),
                          offset = log(PersonYears),
                          data = data)

# Display the summary of the main effects model
print(summary(main_effects_model))

# Check for overdispersion
# Calculate the ratio of residual deviance to degrees of freedom
dispersion_ratio <- main_effects_model$deviance / main_effects_model$df.residual
cat("Dispersion Ratio:", dispersion_ratio, "\n")

if (dispersion_ratio > 1.5) {
  cat("There is evidence of overdispersion (dispersion ratio =", dispersion_ratio, ").\n")
  cat("Consider using a quasipoisson or negative binomial model")
} else {
  cat("No significant overdispersion detected (dispersion ratio =", dispersion_ratio, ").\n")
}

# Adding interaction terms one at a time and compare models using likelihood ratio tests

interaction_models <- list()

model_age_year <- update(main_effects_model, . ~ . + AgeFirstEmployment:YearFirstEmployment)
interaction_models[["Age:Year"]] <- model_age_year

model_age_time <- update(main_effects_model, . ~ . + AgeFirstEmployment:TimeSinceFirstEmployment)
interaction_models[["Age:Time"]] <- model_age_time

#model_year_time <- update(main_effects_model, . ~ . + YearFirstEmployment:TimeSinceFirstEmployment)
#interaction_models[["Year:Time"]] <- model_year_time
#This model was not tested as it was giving this error:
#Error in glm.fit(x = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  : 
#  NA/NaN/Inf in 'x'

anova_results <- lapply(interaction_models, function(model) {
  anova(main_effects_model, model, test = "Chisq")
})

print(anova_results)

# Evaluate p-values to determine if interactions are significant
for (interaction in names(anova_results)) {
  p_value <- anova_results[[interaction]]$`Pr(>Chi)`[2]
  cat("Interaction", interaction, "p-value:", p_value, "\n")
}

# Using the standard p-value threshold of 0.05, we update the model as:
# No interactions should
final_model <- main_effects_model

# Display the summary of the final model
print(summary(final_model))

# Check for overdispersion in the final model
dispersion_ratio_final <- final_model$deviance / final_model$df.residual
dispersion_ratio_final

if (dispersion_ratio_final > 1.5) {
  cat("There is evidence of overdispersion in the final model (dispersion ratio =", dispersion_ratio_final, ").\n")
} else {
  cat("No significant overdispersion detected in the final model (dispersion ratio =", dispersion_ratio_final, ").\n")
}

# Model diagnostics
# Plotting residuals vs. fitted values
png("Residuals_analysis.png")
par(mfrow = c(2, 2))
plot(final_model)
dev.off()

hnp_result <- hnp(
  final_model,
  how.many.out = TRUE,
  paint.out = TRUE,
  warn = FALSE,
  sim = 1000,
  maxit = 1000,
  main = "Half-Normal Plot of Residuals"
)

png("hnp_plot.png")
plot(hnp_result)
dev.off()

influence_measures <- influence.measures(final_model)
summary(influence_measures)

cooksd <- cooks.distance(final_model)

png("cooks_distance.png", width = 800, height = 600)
plot(cooksd, type = "h", main = "Cook's Distance", ylab = "Cook's Distance")
abline(h = 4 / nrow(data), col = "red", lty = 2)
dev.off()


influential_obs <- which(cooksd > (4 / nrow(data)))
cat("Influential Observations based on Cook's Distance:\n")
print(influential_obs)


library(dplyr)

data <- data %>%
  mutate(LogPersonYears = log(PersonYears),
         FittedLogRate = predict(final_model, type = "link"),
         FittedRate = exp(FittedLogRate))

age_rates <- data %>%
  group_by(AgeFirstEmployment) %>%
  summarize(AverageRate = mean(FittedRate))

print("Estimated Cancer Rates by Age at First Employment:")
print(age_rates)

year_rates <- data %>%
  group_by(YearFirstEmployment) %>%
  summarize(AverageRate = mean(FittedRate))

print("Estimated Cancer Rates by Year of First Employment:")
print(year_rates)

time_rates <- data %>%
  group_by(TimeSinceFirstEmployment) %>%
  summarize(AverageRate = mean(FittedRate))

print("Estimated Cancer Rates by Time Since First Employment:")
print(time_rates)

exp_coefficients <- exp(coef(final_model))
conf_intervals <- exp(confint(final_model))

results_table <- data.frame(
  Coefficient = names(exp_coefficients),
  Estimate = exp_coefficients,
  Lower_CI = conf_intervals[, 1],
  Upper_CI = conf_intervals[, 2]
)

print("Exponential of Coefficients with 95% Confidence Intervals:")
print(results_table)

updated_data <- data[-influential_obs, ]

refit_model <- function(updated_data) {
  glm(
    formula(final_model),
    family = poisson(link = "log"),
    data = updated_data,
    offset = log(updated_data$PersonYears)
  )
}

if (length(influential_obs) > 0) {
  final_model_no_influential <- refit_model(updated_data)
  
  original_betas <- exp(coef(final_model))
  betas_no_influential <- exp(coef(final_model_no_influential))
  
  common_coefficients <- intersect(names(original_betas), names(betas_no_influential))
  
  beta_comparison <- data.frame(
    Coefficient = common_coefficients,
    Original = original_betas[common_coefficients],
    No_Influential = betas_no_influential[common_coefficients]
  )
  
  beta_comparison$Abs_Diff <- abs(beta_comparison$Original - beta_comparison$No_Influential)
  beta_comparison$Rel_Abs_Diff <- beta_comparison$Abs_Diff / abs(beta_comparison$Original)
  
  print(beta_comparison)
  
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
  
  p_values_original <- summary(final_model)$coefficients[common_coefficients, "Pr(>|z|)"]
  p_values_no_influential <- summary(final_model_no_influential)$coefficients[common_coefficients, "Pr(>|z|)"]
  
  p_value_comparison <- data.frame(
    Coefficient = common_coefficients,
    Original = p_values_original,
    No_Influential = p_values_no_influential
  )
  
  p_value_long <- melt(p_value_comparison, id.vars = "Coefficient",
                       variable.name = "Model", value.name = "p_value")
  
  p_value_long$minus_log10_pvalue <- -log10(p_value_long$p_value)
  
  max_log10_pvalue <- max(p_value_long$minus_log10_pvalue[is.finite(p_value_long$minus_log10_pvalue)])
  p_value_long$minus_log10_pvalue[is.infinite(p_value_long$minus_log10_pvalue)] <- max_log10_pvalue + 1

  alpha <- 0.05
  num_tests <- nrow(p_value_comparison)
  bonferroni_threshold <- alpha / num_tests
  minus_log10_bonferroni <- -log10(bonferroni_threshold)
  minus_log10_standard <- -log10(alpha)
  
thresholds <- data.frame(
  Threshold = c("Bonferroni Threshold", "Standard Threshold"),
  minus_log10_threshold = c(minus_log10_bonferroni, minus_log10_standard)
)

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

tryCatch({
  hnp_result_no_influential <- hnp(
    final_model_no_influential,
    how.many.out = TRUE,
    paint.out = TRUE,
    warn = FALSE,
    sim = 1000,
    maxit = 1000,
    main = "Half-Normal Plot of Residuals"
  )
  
  png("hnp_plot_no_influential.png")
  plot(hnp_result_no_influential)
  dev.off()
}, error = function(e) {
  print("Envelope cannot be created, showing a potential issue with the model without influential points")
})
  
} else {
  cat("No influential observations detected based on Cook's Distance.\n")
}

sink()