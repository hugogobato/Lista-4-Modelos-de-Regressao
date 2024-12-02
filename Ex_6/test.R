sink()

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
final_model <- glm(final_formula, data = data, family = binomial(link = "logit"))

# Display the summary of the final model
print(summary(final_model))
