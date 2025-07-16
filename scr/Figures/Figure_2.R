# renv::activate()
rm(list=ls())

load("data/data_preprocess/B-splines-normalisation.rda")
# load("data/data_preprocess/target_variables.rda")
# load("data/data_preprocess/fourier-normalisation.rda")
meta_data <- read.csv("data/raw_data/phenotype.csv", sep = ";")

### FIGURE 1.F
cor.test.permu <- function(X, Y){
  r.obs = r.random = NULL
  
  complete_rows <- complete.cases(Y, X)
  Y_clean <- Y[complete_rows, ]
  X_clean <- X[complete_rows, ]
  
  for(b in 1:999){
    boostrap <- sample(nrow(X_clean), replace = T)
    x.boostrap <- as.numeric(dist(scale(X_clean[boostrap,])))
    y.boostrap <- as.numeric(dist(scale(Y_clean[boostrap,])))
    x.permu <- as.numeric(dist(scale(X_clean[sample(nrow(X_clean)),])))
    r.obs = c(r.obs, cor(x.boostrap, y.boostrap))
    r.random = c(r.random, cor(x.permu, as.numeric(dist(scale(Y_clean)))))
  }
  
  return(list(O = r.obs, Null.distri = r.random))
}

blood <- cor.test.permu(meta_data[,8:17], meta_data[,c(3,18,23)])
mir <- cor.test.permu(df$`mir-serum-mean-5months.csv`, meta_data[,c(3,18,23)])
metal <- cor.test.permu(df$`metallomic-5months.csv`, meta_data[,c(3,18,23)])

# Densités
d_metal <- density(metal$O)
d_mir <- density(mir$O)
d_blood <- density(blood$O)
d_null <- density(blood$Null.distri)

pdf(file = "figures/Figure_2A.pdf", width = 5, height = 4)
par(mar=c(4,4,4,2))
# Plot principal
plot(d_metal, col = "black", lwd=2, ylim = c(0,8), xlim = c(-0.25, 0.6),
     main = "Mantel tests between\nvariables describing steatosis and other data", axes=F,
     xlab = "Distribution of Pearson'r correlation\n(999 bootstraps)", cex.main=1)
axis(1)
axis(2)

# Hachures sous la courbe grise
polygon(d_null$x, d_null$y, col = NA, border = NA, density = 20, angle = 45)

# Courbe grise par-dessus les hachures
lines(d_null, col = "grey", lwd = 2, lty=1)

# Autres courbes
lines(d_mir, col = "black", lwd = 2, lty=2)
lines(d_blood, col = "black", lwd = 2, lty=3)

# Legend
legend("topright",
       legend = c("Lipid store vs serum markers", "Lipid store vs MIR", "Lipid store vs Metallome", "Null distribution"),
       col = c("black", "black", "black", "grey"),
       lty = c(1, 2, 3, 1), cex=0.7,
       lwd = 2,
       box.lty = 0)
dev.off()

# NEXT
mir.1 <- cor.test.permu(meta_data[,8:17], df$`mir-serum-mean-5months.csv`)
metal.1 <- cor.test.permu(meta_data[,8:17], df$`metallomic-5months.csv`)

# Densités
d_mir <- density(mir.1$O)
d_metal <- density(metal.1$O)
d_null <- density(metal.1$Null.distri)

pdf(file = "figures/Figure_2B.pdf", width = 5, height = 4)
par(mar=c(4,4,4,2))
# Plot principal
plot(d_mir, col = "black", lwd=2, ylim = c(0,8), xlim = c(-0.25, 0.6),
     main = "Mantel tests between\n serum markers and other data", axes = F, cex.main = 1,
     xlab = "Distribution of Pearson'r correlation\n(999 bootstraps)")
axis(1)
axis(2)

# Courbe rouge
lines(d_metal, col = "black", lwd = 2, lty=2)

# Hachures sous la courbe grise
polygon(d_null$x, d_null$y, col = NA, border = NA, density = 20, angle = 45)

# Courbe grise par-dessus les hachures
lines(d_null, col = "grey", lwd = 2, lty=1)

legend("topright",
       legend = c("MIR vs serum markers", "Metallome vs serum markers", "Null distribution"),
       col = c("black", "black", "grey"),
       lty = c(1, 2, 1), cex=0.7,
       lwd = 2,
       box.lty = 0)
dev.off()