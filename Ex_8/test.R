sink()
updated_data <- data[-influential_obs, ]
final_model_no_influential <- refit_model(data, influential_obs)
print(final_model_no_influential)

