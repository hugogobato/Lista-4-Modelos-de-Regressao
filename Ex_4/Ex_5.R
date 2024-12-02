library(ggplot2)
library(dplyr)
library(MASS)
library(reshape2)
library(hnp)

sink("output.txt")

raw_data <- "
11 low 21 98
11 low 42 96
11 low 62 62
11 medium 21 94
11 medium 42 79
11 medium 62 3
11 high 21 92
11 high 42 41
11 high 62 1
21 low 21 94
21 low 42 93
21 low 62 65
21 medium 21 94
21 medium 42 71
21 medium 62 2
21 high 21 91
21 high 42 30
21 high 62 1
"

data_lines <- unlist(strsplit(raw_data, "\n"))
data_lines <- data_lines[data_lines != ""]

Germ_Temp <- numeric()
Humidity_Level <- character()
Temperature <- numeric()
Germinated <- numeric()

for (line in data_lines) {
  elements <- unlist(strsplit(line, "\\s+"))
  elements <- elements[elements != ""]
  
  Germ_Temp <- c(Germ_Temp, as.numeric(elements[1]))
  Humidity_Level <- c(Humidity_Level, elements[2])
  Temperature <- c(Temperature, as.numeric(elements[3]))
  Germinated <- c(Germinated, as.numeric(elements[4]))
}

data <- data.frame(
  Germ_Temp = factor(Germ_Temp),
  Humidity_Level = factor(Humidity_Level, levels = c("low", "medium", "high")),
  Temperature = factor(Temperature),
  Germinated = Germinated,
  Not_Germinated = 100 - Germinated
)

print(head(data))

print(summary(data))

png("barplot_Germ_Temp.png", width = 800, height = 600)
print(
  ggplot(data, aes(x = Germ_Temp)) +
    geom_bar(fill = "skyblue") +
    labs(title = "Frequency of Germination Temperatures",
         x = "Germination Temperature (°C)",
         y = "Count") +
    theme_bw()
)
dev.off()

png("barplot_Humidity_Level.png", width = 800, height = 600)
print(
  ggplot(data, aes(x = Humidity_Level)) +
    geom_bar(fill = "lightgreen") +
    labs(title = "Frequency of Humidity Levels",
         x = "Humidity Level",
         y = "Count") +
    theme_bw()
)
dev.off()

png("boxplot_Germinated_by_Germ_Temp.png", width = 800, height = 600)
print(
  ggplot(data, aes(x = Germ_Temp, y = Germinated)) +
    geom_boxplot(fill = "orange") +
    labs(title = "Germinated Seeds by Germination Temperature",
         x = "Germination Temperature (°C)",
         y = "Number of Germinated Seeds") +
    theme_bw()
)
dev.off()

png("boxplot_Germinated_by_Humidity_Level.png", width = 800, height = 600)
print(
  ggplot(data, aes(x = Humidity_Level, y = Germinated)) +
    geom_boxplot(fill = "violet") +
    labs(title = "Germinated Seeds by Humidity Level",
         x = "Humidity Level",
         y = "Number of Germinated Seeds") +
    theme_bw()
)
dev.off()

png("interaction_plot_Germinated_Temperature_Humidity.png", width = 800, height = 600)
print(
  ggplot(data, aes(x = Temperature, y = Germinated, color = Humidity_Level, group = Humidity_Level)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    facet_wrap(~ Germ_Temp) +
    labs(title = "Interaction of Temperature and Humidity Level on Germination",
         x = "Temperature Level (°C)",
         y = "Number of Germinated Seeds",
         color = "Humidity Level") +
    theme_bw()
)
dev.off()

png("interaction_plot_Germinated_Germ_Temp_Humidity.png", width = 800, height = 600)
print(
  ggplot(data, aes(x = Germ_Temp, y = Germinated, color = Humidity_Level, group = Humidity_Level)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    facet_wrap(~ Temperature) +
    labs(title = "Interaction of Germination Temperature and Humidity Level on Germination",
         x = "Germination Temperature (°C)",
         y = "Number of Germinated Seeds",
         color = "Humidity Level") +
    theme_bw()
)
dev.off()


full_model <- glm(
  cbind(Germinated, Not_Germinated) ~ (Germ_Temp + Humidity_Level + Temperature)^2,
  family = binomial(link = "logit"),
  data = data
)

selected_model <- stepAIC(
  full_model,
  direction = "both",
  trace = FALSE
)

print(summary(selected_model))

coefficients <- coef(selected_model)
odds_ratios <- exp(coefficients)
cat("Odds Ratios of the coefficients:")
print(odds_ratios)

residual_deviance <- deviance(selected_model)
df_residual <- df.residual(selected_model)
overdispersion_ratio <- residual_deviance / df_residual
cat("Overdispersion Ratio:", overdispersion_ratio, "\n")

# Since the data are binomial counts, overdispersion can occur
# If overdispersion_ratio > 1, consider using quasibinomial family


hnp_result <- hnp(
  selected_model,
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


cooksd <- cooks.distance(selected_model)

png("cooks_distance.png", width = 800, height = 600)
plot(cooksd, type = "h", main = "Cook's Distance", ylab = "Cook's Distance")
abline(h = 4 / nrow(data), col = "red", lty = 2)
dev.off()


influential_obs <- which(cooksd > (4 / nrow(data)))
cat("Influential Observations based on Cook's Distance:\n")
print(influential_obs)

refit_model <- function(data, exclude_rows) {
  updated_data <- data[-exclude_rows, ]
  glm(
    formula(selected_model),
    family = binomial(link = "logit"),
    data = updated_data
  )
}

if (length(influential_obs) > 0) {
  selected_model_no_influential <- refit_model(data, influential_obs)
  
  original_betas <- coef(selected_model)
  betas_no_influential <- coef(selected_model_no_influential)
  
  common_coefficients <- c(
  "(Intercept)",
  "Germ_Temp21",
  "Humidity_Levelmedium",
  "Humidity_Levelhigh",
  "Temperature42",
  "Temperature62",
  "Humidity_Levelmedium:Temperature42",
  "Humidity_Levelhigh:Temperature42",
  "Humidity_Levelmedium:Temperature62"
)
  
  original_betas_common <- original_betas[common_coefficients]
  betas_no_influential_common <- betas_no_influential[common_coefficients]
  

  beta_comparison <- data.frame(
    Coefficient = common_coefficients,
    Original = original_betas_common,
    No_Influential = betas_no_influential_common
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
  

  p_values_original <- summary(selected_model)$coefficients[common_coefficients, "Pr(>|z|)"]
  p_values_no_influential <- summary(selected_model_no_influential)$coefficients[common_coefficients, "Pr(>|z|)"]
  
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
  minus_log10_threshold <- -log10(bonferroni_threshold)
  
  png("p_values_comparison_plot.png", width = 1000, height = 600)
  p_pvalues <- ggplot(p_value_long, aes(x = Coefficient, y = minus_log10_pvalue, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    geom_hline(yintercept = minus_log10_threshold, linetype = "dashed", color = "red", size = 1) +
    labs(
      x = "Coefficient",
      y = expression(-log[10](p~value)),
      title = "Comparison of p-values Across Models (Common Coefficients)",
      fill = "Model"
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    ) +
    annotate(
      "text",
      x = 1,
      y = minus_log10_threshold + 0.5,
      label = paste("Bonferroni Threshold (p =", signif(bonferroni_threshold, 3), ")"),
      color = "red",
      hjust = 0,
      size = 5
    )
  print(p_pvalues)
  dev.off()
  
  print(p_value_comparison)
  
} else {
  cat("No influential observations detected based on Cook's Distance.\n")
}
sink()