sink()

data$dose <- factor(data$dose, levels = c(0, 80, 160, 235, 310))

ggplot(data, aes(x = dose, y = tovos)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  labs(
    title = "Total Number of Hatched Eggs by Nitrofen Dose",
    x = "Nitrofen Dose (mg/L)",
    y = "Total Number of Hatched Eggs"
  ) +
  theme_classic()

# Save the boxplot
ggsave("boxplot_tovos_vs_dose.png")
sink()