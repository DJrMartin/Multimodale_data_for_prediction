# renv::activate()
rm(list=ls())
library(FactoMineR)

load("/data/data_preprocess/B-splines-normalisation.rda")
# load("data/data_preprocess/target_variables.rda")
# load("data/data_preprocess/fourier-normalisation.rda")
meta_data <- read.csv("data/raw_data/phenotype.csv", sep = ";")

col.group = as.character(factor(meta_data$Group, levels(as.factor(meta_data$Group)), c("#AFB9C8", "#8497B0", "#F4B183", "#F1CDB1")))
pch.group = as.numeric(as.character(factor(meta_data$Group, levels(as.factor(meta_data$Group)), c(15:18))))

# M vs C
par(mar=c(4,4,4,2))
plot(meta_data$Fsp27~meta_data$Triglycerides.Hp, pch = pch.group,  col = col.group, main = "Relationship between cellular and tissular\nindicators of Steatosis",
     axes = FALSE, xlim=c(0, 100), xlab = "Hepatic Trigly. (mg/g)", ylab = "Histological Steatosis (%)", cex.main=1) 
axis(1)
axis(2, at=seq(0,1,0.2))
legend("topleft", legend = c("CTL", "HFHC", "IRON", "IRON+HFHC"), col = c("#AFB9C8", "#8497B0","#F1CDB1","#F4B183"), ncol=2, 
       bty='n', cex = 0.8, pch = c(15, 16, 18, 17))

# M vs T
par(mar=c(4,4,4,2))
plot(meta_data$Steatosis~meta_data$Triglycerides.Hp, pch = pch.group,  col = col.group, main = "Relationship between cellular and tissular\nindicators of Steatosis",
     axes = FALSE, ylim=c(0,1.2), xlim=c(0, 100), xlab = "Hepatic Trigly. (mg/g)", ylab = "Histological Steatosis (%)", cex.main=1) 
axis(1)
axis(2, at=seq(0,1,0.2))
legend("topleft", legend = c("CTL", "HFHC", "IRON", "IRON+HFHC"), col = c("#AFB9C8", "#8497B0","#F1CDB1","#F4B183"), ncol=2, 
       bty='n', cex = 0.8, pch = c(15, 16, 18, 17))

group = factor(meta_data$Group, c("CON", "HF", "IRON", "HF+IRON"))
boxplot(meta_data$Steatosis~group, col = c("#AFB9C8", "#8497B0","#F1CDB1","#F4B183"), 
        xlab = "", ylab = "", axes = F)
axis(1)
axis(2)

colnames(meta_data[, c(2,3,5,11,12,18,23)])
res.pca <- PCA(meta_data[, c(2,3,5,18,23)], quali.sup = c(1,3), graph=T)

plot(res.pca$ind$coord, col = col.group, pch = pch.fibrosis, main = "PCA on variables describing steatosis",  cex.main=1)
abline(v=0, lty="dashed", col='grey')
abline(h=0, lty="dashed", col='grey')
legend("topright", legend = c("CTL", "HFHC", "IRON", "IRON+HFHC"),fill= c("#D6DCE5", "#8497B0","#FBE5D6","#F4B183"), ncol=2, 
       bty='n', cex=0.8)
legend("bottomright", legend = c("0", "1"), pch = c(16, 17),col = c("grey"), 
       bty='n', cex=0.8, title="Fibrosis score")

DF = scale(meta_data[,c(9,10,13,16,17, 
                            3, 18, 23, 
                            5, 20, 
                            21, 22, 24,
                            14,15, 19)])
n = nrow(DF)
p = ncol(DF)
  
B = 1000
PC1_exp = PC2_exp = NULL
for (i in 1:p){
  coef_lm_PC1=coef_lm_PC2=NULL
  for (cnt in 1:B){
    boostrap=sample(1:50, 50, replace = T)
    coef_lm_PC1 = c(coef_lm_PC1,lm(scale(res.pca$ind$coord[boostrap,1])~scale(DF[boostrap,i]))$coefficients[-1])
    coef_lm_PC2 = c(coef_lm_PC2,lm(scale(res.pca$ind$coord[boostrap,2])~scale(DF[boostrap,i]))$coefficients[-1])
  }
  PC1_exp = cbind(PC1_exp,coef_lm_PC1)
  PC2_exp = cbind(PC2_exp,coef_lm_PC2)
}

par(mar=c(7,5,4,1))
plot(NA, xlim = c(0.5,(p+.5)), ylim = c(-1,1), main = "Explication of the dimensions of the PCA in (C)", cex.main = 1,
     axes = F, xlab = "", ylab = expression(paste("95% confidence interval of\nthe Pearson's r correlation ")))
abline(h=0, col="black", lty='dashed')
axis(2)
## DIM1
col.var = rep("grey", 16)
col.var[c(2, 4, 6:8)] = "#8497B0"

bornes_inf = apply(PC1_exp,2,quantile,probs=.025)
bornes_sup = apply(PC1_exp,2,quantile,probs=.975)
medianes = apply(PC1_exp,2,quantile,probs=.5)
X.Y=data.frame(x1=0.8:(p-.2), y1=bornes_inf,
               x2=0.8:(p-.2), y2=bornes_sup)
X.Y=X.Y[which(is.na(X.Y[,2])==F),]
segments(X.Y$x1, X.Y$y1, X.Y$x2, X.Y$y2, lwd=2)
points(x=0.8:(p-.2),y=medianes, lwd=2, cex=1.5 , pch=15, col=col.var)
## DIM2
col.var = rep("grey", 16)
col.var[c(1, 9, 14:15)] = "#F4B183"
bornes_inf = apply(PC2_exp,2,quantile,probs=.025)
bornes_sup = apply(PC2_exp,2,quantile,probs=.975)
medianes = apply(PC2_exp,2,quantile,probs=.5)
X.Y=data.frame(x1=1.2:(p+.2), y1=bornes_inf,
               x2=1.2:(p+.2), y2=bornes_sup)
X.Y=X.Y[which(is.na(X.Y[,2])==F),]
segments(X.Y$x1, X.Y$y1, X.Y$x2, X.Y$y2, lwd=2)
points(x=1.2:(p+.2),y=medianes, lwd=3, cex=1.5 , pch=16, col=col.var)

box()
abline(v=c(5.5,8.5,10.5,13.5), col="gray", lwd=2)
axis(1,1:p,labels=FALSE)
names.var = c("PAL","ALAT","ASAT","ASAT/ALAT" ,"CK","Steatosis (%)","Hepatic Trigly (mg/g)","Fsp27 mRNA",           
              "Fibrosis score","Col1a mRNA","NRF2 mRNA","A2M mRNA","Crp mRNA","Total serum iron","Transferin Sat.","Hamp mRNA")
text(0.8:p,rep(-1.3, p), names.var, 
     srt = 45, xpd=NA, adj=c(1,1), cex = 0.8)
legend("topleft", legend = c("r correlation with Dim.1", "r correlation with Dim.2"),
       bty = 'n', pch=c(15, 16), cex=0.8, title.adj = 0)
text(x=c(1,6,9,11,14), y=-0.95, c("(1)", "(2)", "(3)", "(4)", "(5)"))
text(x=c(7.1), y=-0.95, "Used for\nthe PCA", cex=0.7)

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
       legend = c("Lipid store vs blood markers", "Lipid store vs MIR", "Lipid store vs Metallome", "Null distribution"),
       col = c("black", "black", "black", "grey"),
       lty = c(1, 2, 3, 1), cex=0.7,
       lwd = 2,
       box.lty = 0)

# NEXT
mir.1 <- cor.test.permu(meta_data[,8:17], df$`mir-serum-mean-5months.csv`)
metal.1 <- cor.test.permu(meta_data[,8:17], df$`metallomic-5months.csv`)

# Densités
d_mir <- density(mir.1$O)
d_metal <- density(metal.1$O)
d_null <- density(metal.1$Null.distri)

# Plot principal
plot(d_mir, col = "black", lwd=2, ylim = c(0,8), xlim = c(-0.25, 0.6),
     main = "Mantel tests between\n blood markers and other data", axes = F, cex.main = 1,
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
       legend = c("MIR vs blood markers", "Metallome vs blood markers", "Null distribution"),
       col = c("black", "black", "grey"),
       lty = c(1, 2, 1), cex=0.7,
       lwd = 2,
       box.lty = 0)

