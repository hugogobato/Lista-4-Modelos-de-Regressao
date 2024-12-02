sink()

selected_model <- glm(final_formula, data = data, family = binomial(link = "logit"))

print(summary(selected_model))

print(final_formula)