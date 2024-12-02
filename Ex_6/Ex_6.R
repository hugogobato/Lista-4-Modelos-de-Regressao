library(MASS)        
library(hnp)         
library(ggplot2)     
library(reshape2)    

data_url <- "https://www.ime.usp.br/~giapaula/leuce.txt"
data <- read.table(data_url, header = FALSE)

sink("output.txt")

colnames(data) <- c("age", "differential_count", "bone_marrow_infiltration", 
                    "leukemia_cells", "malignancy", "max_temp_before_treatment",
                    "treatment_efficiency", "survival_time", "status")

data$treatment_efficiency <- as.numeric(as.character(data$treatment_efficiency))

print(head(data))

explanatory_vars <- c("age", "differential_count", "bone_marrow_infiltration", 
                      "leukemia_cells", "malignancy", "max_temp_before_treatment")

interaction_terms <- combn(explanatory_vars, 2, FUN = function(x) paste(x, collapse = ":"))

# Set thresholds for entry (PE) and stay (PS)
entry_threshold <- 0.20
stay_threshold <- 0.20

# Define base formula (intercept-only model)
base_formula <- as.formula("treatment_efficiency ~ 1")

# Extract explanatory variables and interaction terms
all_predictors <- c(explanatory_vars, interaction_terms)

# Initialize variables for tracking
current_formula <- base_formula
remaining_predictors <- all_predictors
current_model <- glm(current_formula, data = data, family = binomial(link = "logit"))

# Forward stepwise selection
step <- 1
selected_terms <- c()  # Keep track of selected terms
repeat {
  # Fit models with one additional predictor at a time
  candidate_models <- lapply(remaining_predictors, function(pred) {
    formula_string <- paste(deparse(current_formula), "+", pred)
    formula <- as.formula(formula_string)
    glm(formula, data = data, family = binomial(link = "logit"))
  })
  
  # Extract p-values for the new terms in each candidate model
  p_values <- sapply(seq_along(candidate_models), function(i) {
    model <- candidate_models[[i]]
    coef_summary <- summary(model)$coefficients
    added_terms <- setdiff(rownames(coef_summary), rownames(summary(current_model)$coefficients))
    if (length(added_terms) > 0) {
      coef_summary[added_terms, "Pr(>|z|)"]
    } else {
      NA
    }
  })
  
  # Assign predictor names to the p_values vector
  names(p_values) <- remaining_predictors
  
  # Identify the best predictor (lowest p-value within threshold)
  best_p_value <- min(p_values, na.rm = TRUE)
  if (best_p_value <= entry_threshold) {
    best_predictor <- names(which.min(p_values))
    current_formula <- as.formula(paste(deparse(current_formula), "+", best_predictor))
    current_model <- glm(current_formula, data = data, family = binomial(link = "logit"))
    remaining_predictors <- setdiff(remaining_predictors, best_predictor)
    selected_terms <- c(selected_terms, best_predictor)
    cat(sprintf("Step %d: Added '%s' with p-value = %.4f\n", step, best_predictor, best_p_value))
    step <- step + 1
  } else {
    cat("No additional predictors meet the entry threshold. Stopping forward selection.\n")
    break
  }
}

# Backward elimination for terms with p-value > stay threshold
final_model <- current_model
repeat {
  coef_summary <- summary(final_model)$coefficients
  p_values <- coef_summary[-1, "Pr(>|z|)"]  # Exclude intercept
  
  if (length(p_values) > 0 && max(p_values) > stay_threshold) {
    # Identify term with highest p-value
    worst_term <- names(which.max(p_values))
    worst_p_value <- max(p_values)
    
    # Remove the term only if it is not part of an interaction or required for another term
    interaction_terms <- grep(":", attr(terms(final_model), "term.labels"), value = TRUE)
    associated_terms <- unique(unlist(strsplit(interaction_terms, ":")))
    if (!(worst_term %in% associated_terms)) {
      updated_terms <- setdiff(attr(terms(final_model), "term.labels"), worst_term)
      current_formula <- as.formula(paste("treatment_efficiency ~", paste(updated_terms, collapse = " + ")))
      final_model <- glm(current_formula, data = data, family = binomial(link = "logit"))
      selected_terms <- setdiff(selected_terms, worst_term)
      cat(sprintf("Removed '%s' with p-value = %.4f\n", worst_term, worst_p_value))
    } else {
      cat(sprintf("Skipped removing '%s' as it is part of an interaction.\n", worst_term))
      break
    }
  } else {
    cat("All remaining terms meet the stay threshold. Final model selected.\n")
    break
  }
}

# Combine the selected terms into the final formula
final_formula_string <- paste("treatment_efficiency ~", paste(selected_terms, collapse = " + "))

# Convert the string to a formula
final_formula <- as.formula(final_formula_string)

# Fit the final model
selected_model <- glm(final_formula, data = data, family = binomial(link = "logit"))

print(summary(selected_model))

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
  
  common_coefficients <- intersect(names(original_betas), names(betas_no_influential))
  
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

exp_coef <- exp(coef(selected_model))
exp_confint <- exp(confint(selected_model))

odds_ratios <- data.frame(
  Coefficient = names(coef(selected_model)),
  Odds_Ratio = exp_coef,
  Lower_CI = exp_confint[,1],
  Upper_CI = exp_confint[,2]
)

print("Odds Ratios and 95% Confidence Intervals:")
print(odds_ratios)

# Interpretation example:
# For a continuous variable like "age", an odds ratio greater than 1 indicates that as age increases,
# the odds of treatment being satisfactory increase. An odds ratio less than 1 indicates the opposite.
