library(ggplot2)
library(dplyr)
library(MASS)


url <- "https://www.ime.usp.br/~giapaula/claims.txt"
data <- read.table(url, header = FALSE)

colnames(data) <- c("valorv", "expos", "nsinistros", "csinistros", "tipov",
                    "idadev", "sexoc", "areac", "idadec")


png("histograms_numerical.png")
par(mfrow = c(2, 2))
hist(data$valorv, main = "Vehicle Value", xlab = "Value (in 10,000 AUD)")
hist(data$expos, main = "Exposure", xlab = "Exposure")
hist(data$nsinistros, main = "Number of Claims", xlab = "Number of Claims")
hist(data$csinistros, main = "Total Claim Cost", xlab = "Cost (AUD)")
dev.off()

data <- data %>%
  mutate(cmsinistros = ifelse(nsinistros > 0, csinistros / nsinistros, NA))

data_clean <- data %>% filter(!is.na(cmsinistros))

print(summary(data[, c("valorv", "expos", "cmsinistros")]))

png("histogram_cmsinistros.png")
hist(data_clean$cmsinistros, main = "Average Cost per Claim",
     xlab = "Cost per Claim (AUD)")
dev.off()

boxplot_plot <- ggplot(data_clean, aes(x = factor(idadev), y = cmsinistros)) +
  geom_boxplot() +
  xlab("Vehicle Age Category") +
  ylab("Cmsinistros")

ggsave("boxplot_idadev.png", plot = boxplot_plot, width = 8, height = 6, dpi = 300)

boxplot_plot <- ggplot(data_clean, aes(x = factor(tipov), y = cmsinistros)) +
  geom_boxplot() +
  xlab("Car Type Category") +
  ylab("Cmsinistros")

ggsave("boxplot_tipov.png", plot = boxplot_plot, width = 8, height = 6, dpi = 300)

boxplot_plot <- ggplot(data_clean, aes(x = factor(sexoc), y = cmsinistros)) +
  geom_boxplot() +
  xlab("Gender") +
  ylab("Cmsinistros")

ggsave("boxplot_sexoc.png", plot = boxplot_plot, width = 8, height = 6, dpi = 300)

boxplot_plot <- ggplot(data_clean, aes(x = factor(areac), y = cmsinistros)) +
  geom_boxplot() +
  xlab("Living Area Category") +
  ylab("Cmsinistros")

ggsave("boxplot_areac.png", plot = boxplot_plot, width = 8, height = 6, dpi = 300)

boxplot_plot <- ggplot(data_clean, aes(x = factor(idadec), y = cmsinistros)) +
  geom_boxplot() +
  xlab("Driver's Age Category") +
  ylab("Cmsinistros")

ggsave("boxplot_idadec.png", plot = boxplot_plot, width = 8, height = 6, dpi = 300)

scatter1 <- ggplot(data_clean, aes(x = valorv, y = cmsinistros)) +
  geom_point(alpha = 0.6) +
  labs(title = "Scatter Plot: cmsinistros ~ valorv",
       x = "Vehicle Value (valorv)",
       y = "Average Claim Cost (cmsinistros)") +
  theme_classic()

ggsave("scatter_cmsinistros_valorv.png", plot = scatter1, width = 8, height = 6, dpi = 300)

scatter2 <- ggplot(data_clean, aes(x = valorv, y = log(cmsinistros))) +
  geom_point(alpha = 0.6) +
  labs(title = "Scatter Plot: log(cmsinistros) ~ valorv",
       x = "Vehicle Value (valorv)",
       y = "Log of Average Claim Cost (log(cmsinistros))") +
  theme_classic()

ggsave("scatter_log_cmsinistros_valorv.png", plot = scatter2, width = 8, height = 6, dpi = 300)

scatter3 <- ggplot(data_clean, aes(x = expos, y = cmsinistros)) +
  geom_point(alpha = 0.6) +
  labs(title = "Scatter Plot: cmsinistros ~ expos",
       x = "Exposure (expos)",
       y = "Average Claim Cost (cmsinistros)") +
  theme_classic()

ggsave("scatter_cmsinistros_expos.png", plot = scatter3, width = 8, height = 6, dpi = 300)

scatter4 <- ggplot(data_clean, aes(x = expos, y = log(cmsinistros))) +
  geom_point(alpha = 0.6) +
  labs(title = "Scatter Plot: log(cmsinistros) ~ expos",
       x = "Exposure (expos)",
       y = "Log of Average Claim Cost (log(cmsinistros))") +
  theme_classic()

ggsave("scatter_log_cmsinistros_expos.png", plot = scatter4, width = 8, height = 6, dpi = 300)


gamma_model <- glm(cmsinistros ~ valorv + expos + factor(tipov) + factor(idadev) +
                     factor(sexoc) + factor(areac) + factor(idadec),
                   family = Gamma(link = "log"), data = data_clean)

# Inverse Gaussian model is not converging. I am getting the error:
# Error: inner loop 1; cannot correct step size
#In addition: Warning messages:
#1: step size truncated due to divergence
#2: step size truncated due to divergence   

# Possible explanations:
# 1. Extreme Variance of cmsinistros
# 2. Outliers or Skewness of cmsinistros: Inverse Gaussian models are sensitive to such skewness and outliers, making convergence difficult.
# 3. Small Exposure Values (expos)
# 4. Convergence Sensitivity: The Inverse Gaussian likelihood maximization is highly sensitive to starting values, small changes in predictors, or high correlations among predictors.

# So I will: 1.Address Outliers, 2. Rescale Predictors, 3. Adjust the Iteration Parameters

data_clean <- data_clean %>% 
  filter(cmsinistros < quantile(cmsinistros, 0.99))  # Remove top 1% of values

data_clean <- data_clean %>%
  mutate(valorv_scaled = scale(valorv), expos_scaled = scale(expos))


invgauss_model <- glm(cmsinistros ~ valorv + expos + factor(tipov) + factor(idadev) +
                     factor(sexoc) + factor(areac) + factor(idadec),
                      family = inverse.gaussian(link = "log"), data = data_clean,control = glm.control(maxit = 100, epsilon = 1e-8))

# It worked

png("gamma_model_diagnostics.png")
par(mfrow = c(2, 2))
plot(gamma_model)
dev.off()

png("invgauss_model_diagnostics.png")
par(mfrow = c(2, 2))
plot(invgauss_model)
dev.off()

png("gamma_residuals_vs_fitted.png")
plot(gamma_model$fitted.values, residuals(gamma_model, type = "deviance"),
     main = "Gamma Model Residuals vs Fitted",
     xlab = "Fitted Values", ylab = "Deviance Residuals")
abline(h = 0, col = "red")
dev.off()

png("invgauss_residuals_vs_fitted.png")
plot(invgauss_model$fitted.values, residuals(invgauss_model, type = "deviance"),
     main = "Inverse Gaussian Model Residuals vs Fitted",
     xlab = "Fitted Values", ylab = "Deviance Residuals")
abline(h = 0, col = "red")
dev.off()

gamma_dispersion <- sum(residuals(gamma_model, type = "pearson")^2) /
  gamma_model$df.residual
invgauss_dispersion <- sum(residuals(invgauss_model, type = "pearson")^2) /
  invgauss_model$df.residual

writeLines(c(
  paste("Gamma Model Dispersion:", gamma_dispersion),
  paste("Inverse Gaussian Model Dispersion:", invgauss_dispersion)
), con = "model_dispersion.txt")

output_dir <-"C:/Users/hugog/Desktop/Master_Courses/Regression_Models/List 4"

mu_hats_gamma <- predict(gamma_model, type = "response")
mu_hats_invgauss <- predict(invgauss_model, type = "response")

plot_data <- data.frame(Index = 1:length(mu_hats_gamma), Predictions = mu_hats_gamma )

p <- ggplot(plot_data, aes(x = Index, y = Predictions)) +
  geom_point(color = "blue") +
  labs(title = "Scatter Plot of Predicted Values (mu_hats)",
       x = "Observation Index",
       y = "Predicted Values (mu_hats)") +
   theme_bw()


ggsave(filename = "gamma_mu_hats_scatter_plot.png", plot = p)

plot_data <- data.frame(Index = 1:length(mu_hats_invgauss), Predictions = mu_hats_invgauss)

p <- ggplot(plot_data, aes(x = Index, y = Predictions)) +
  geom_point(color = "blue") +
  labs(title = "Scatter Plot of Predicted Values (mu_hats)",
       x = "Observation Index",
       y = "Predicted Values (mu_hats)") +
   theme_bw()


ggsave(filename = "invgauss_mu_hats_scatter_plot.png", plot = p)

X <- model.matrix(gamma_model)
w <- gamma_model$weights
W <- diag(w)
fitted_values <- fitted(gamma_model)

# Leverage Points
auxh <- solve(t(X) %*% W %*% X)
H <- sqrt(W) %*% X %*% auxh %*% t(X) %*% sqrt(W)
h <- diag(H)


png(filename = file.path(output_dir, "gamma_leverage_points.png"))
plot(fitted_values, h, xlab = "Fitted Values", ylab = "Leverage", main = "Leverage Points")
identify(fitted_values, h, n = 3)
dev.off()

# To avoid dividing by zero
h <- ifelse(h >= 1, 0.975, h)
# Deviance Residuals
aux_tdi <- resid(gamma_model, type = 'deviance')
fi <- 1 / summary(gamma_model)$dispersion
tdi <- aux_tdi * sqrt(fi / (1 - h))

png(filename = file.path(output_dir, "gamma_deviance_residuals.png"))
plot(fitted_values, tdi, xlab = "Fitted Values", ylab = "Deviance Residuals", main = "Deviance Residuals")
identify(fitted_values, tdi, n = 3)
dev.off()

# Pearson Residuals
aux_tsi <- resid(gamma_model, type = 'pearson')
tsi <- aux_tsi * sqrt(fi / (1 - h))

png(filename = file.path(output_dir, "gamma_pearson_residuals.png"))
plot(fitted_values, tsi, xlab = "Fitted Values", ylab = "Pearson Residuals", main = "Pearson Residuals")
identify(fitted_values, tsi, n = 2)
dev.off()

# Cook's Distance
ldi <- h * (tsi^2) / (1 - h)

png(filename = file.path(output_dir, "gamma_cooks_distance.png"))
plot(ldi, xlab = "Index", ylab = "Cook's Distance", main = "Cook's Distance")
identify(ldi, n = 2)
dev.off()

sink("Output_Model_Diagnostics_Ex_2.txt")

cat("Gamma Model:", "\n")
print(summary(gamma_model))
print(exp(coefficients(gamma_model)))
print(
list(
  Leverage = h,
  Deviance_Residuals = tdi,
  Pearson_Residuals = tsi,
  Cooks_Distance = ldi
))


X <- model.matrix(invgauss_model)
w <- invgauss_model$weights
W <- diag(w)
fitted_values <- fitted(invgauss_model)

# Leverage Points
auxh <- solve(t(X) %*% W %*% X)
H <- sqrt(W) %*% X %*% auxh %*% t(X) %*% sqrt(W)
h <- diag(H)


png(filename = file.path(output_dir, "invgauss_leverage_points.png"))
plot(fitted_values, h, xlab = "Fitted Values", ylab = "Leverage", main = "Leverage Points")
identify(fitted_values, h, n = 3)
dev.off()

# To avoid dividing by zero
h <- ifelse(h >= 1, 0.975, h)
# Deviance Residuals
aux_tdi <- resid(invgauss_model, type = 'deviance')
fi <- 1 / summary(invgauss_model)$dispersion
tdi <- aux_tdi * sqrt(fi / (1 - h))

png(filename = file.path(output_dir, "invgauss_deviance_residuals.png"))
plot(fitted_values, tdi, xlab = "Fitted Values", ylab = "Deviance Residuals", main = "Deviance Residuals")
identify(fitted_values, tdi, n = 3)
dev.off()

# Pearson Residuals
aux_tsi <- resid(invgauss_model, type = 'pearson')
tsi <- aux_tsi * sqrt(fi / (1 - h))

png(filename = file.path(output_dir, "invgauss_pearson_residuals.png"))
plot(fitted_values, tsi, xlab = "Fitted Values", ylab = "Pearson Residuals", main = "Pearson Residuals")
identify(fitted_values, tsi, n = 2)
dev.off()

# Cook's Distance
ldi <- h * (tsi^2) / (1 - h)

png(filename = file.path(output_dir, "invgauss_cooks_distance.png"))
plot(ldi, xlab = "Index", ylab = "Cook's Distance", main = "Cook's Distance")
identify(ldi, n = 2)
dev.off()

cat("Inverse Gaussian Model:", "\n")
print(summary(invgauss_model))
print(exp(coefficients(invgauss_model)))
print(
list(
  Leverage = h,
  Deviance_Residuals = tdi,
  Pearson_Residuals = tsi,
  Cooks_Distance = ldi
))