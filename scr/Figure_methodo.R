# library
library(ggplot2)

# sigmoïde
sigmoid <- function(x, x0 = 0, k = 1) {
  1 / (1 + exp(-k * (x - x0)))
}

# data
x <- seq(-10, 30, length.out = 500)
df <- data.frame(
  x = rep(x, 3),
  y = c(sigmoid(x, x0 = 0, k = 1),
        sigmoid(x, x0 = 10, k = 1),
        sigmoid(x, x0 = 20, k = 1)),
  scale = factor(rep(c("Molecular", "Cellular", "Tissular"), each = length(x)),
                 levels = c("Molecular", "Cellular", "Tissular"))
)

line_types <- c(1:3)

line_colors <- c("Molecular" = "black", 
                 "Cellular" = "black", 
                 "Tissular" = "black")

# figures
pdf(file = "figures/methods_A.pdf", width = 5.5, height = 4)
ggplot(df, aes(x = x, y = y, color = scale, linetype = scale)) +
  geom_line(size = 1) +
  scale_color_manual(values = line_colors) +
  scale_linetype_manual(values = line_types) +
  labs(
    x = "Conceptual axes (e.g., time, scale, severity)",
    y = "Activation / Response",
    color = "Scale",
    linetype = "Scale"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top") +
  # Ajout des textes sur les courbes
  annotate("text", x = -3, y = 0.20, label = "mRNA level of Fsp27", color = "black", angle = 78, hjust = 0) +
  annotate("text", x = 7, y = 0.20, label = "Hepatic triglycerides", color = "black", angle = 78, hjust = 0) +
  annotate("text", x = 17, y = 0.20, label = "Histological steatosis", color = "black", angle = 78, hjust = 0)
dev.off()
