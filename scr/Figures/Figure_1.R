# renv::activate()
rm(list=ls())
library(FactoMineR)

load("data/data_preprocess/B-splines-normalisation.rda")
# load("data/data_preprocess/target_variables.rda")
# load("data/data_preprocess/fourier-normalisation.rda")
meta_data <- read.csv("data/raw_data/phenotype.csv", sep = ";")

col.group = as.character(factor(meta_data$Group, levels(as.factor(meta_data$Group)), c("#AFB9C8", "#8497B0", "#F4B183", "#F1CDB1")))
pch.group = as.numeric(as.character(factor(meta_data$Group, levels(as.factor(meta_data$Group)), c(15:18))))

################# FIGURE 1A-C ######################
pdf(file = "figures/Figures_1A-C.pdf", width = 8, height = 4)
layout(matrix(c(1:3), nrow = 1))
# Molecular
group = factor(meta_data$Group, c("CON", "HF", "IRON", "HF+IRON"))
boxplot(meta_data$Fsp27~group, col = c("#AFB9C8", "#8497B0","#F1CDB1","#F4B183"), main = 'Molecular scale',
        xlab = "", ylab = "mRNA level of Fsp27 (CT)", axes = F, ylim = c(-10,3))
legend("topleft", legend = c("CTL", "HFHC", "IRON", "IRON+HFHC"), fill = c("#AFB9C8", "#8497B0","#F1CDB1","#F4B183"), ncol=2, 
       bty='n', cex = 0.8)
axis(2)
text(c(2,4), y = c(1,1), c("a", "b"))

# Cellular
group = factor(meta_data$Group, c("CON", "HF", "IRON", "HF+IRON"))
boxplot(meta_data$Triglycerides.Hp~group, col = c("#AFB9C8", "#8497B0","#F1CDB1","#F4B183"), main = 'Cellular scale',
        xlab = "", ylab = "Hepatic triglycerides (mg/g)", axes = F, ylim = c(0,120))
legend("topleft", legend = c("CTL", "HFHC", "IRON", "IRON+HFHC"), fill = c("#AFB9C8", "#8497B0","#F1CDB1","#F4B183"), ncol=2, 
       bty='n', cex = 0.8)
axis(2)
text(c(2,4), y = c(100), c("a", "b"))

# Tissular
group = factor(meta_data$Group, c("CON", "HF", "IRON", "HF+IRON"))
boxplot(meta_data$Steatosis~group, col = c("#AFB9C8", "#8497B0","#F1CDB1","#F4B183"), main = 'Tissular scale',
        xlab = "", ylab = "Hepatic histological steatosis (%)", axes = F, ylim = c(0,1.2))
legend("topleft", legend = c("CTL", "HFHC", "IRON", "IRON+HFHC"), fill = c("#AFB9C8", "#8497B0","#F1CDB1","#F4B183"), ncol=2, 
       bty='n', cex = 0.8)
axis(2, at = seq(0, 1, by = 0.2))
text(c(2,4), y = c(1.05), c("a", "b"))
dev.off()

################# FIGURE 1D ######################
pdf(file = "figures/Figures_1D.pdf", width = 5, height = 4)
par(mar=c(4,4,4,2))
plot(meta_data$Steatosis~meta_data$Triglycerides.Hp, pch = pch.group,  col = col.group, main = "Relationship between cellular and tissular\nindicators of Steatosis",
     axes = FALSE, ylim=c(0,1.3), xlim=c(0, 100), xlab = "Hepatic Trigly. (mg/g)", ylab = "Histological Steatosis (%)", cex.main=1) 
axis(1)
axis(2, at=seq(0, 1, 0.2))
legend("topright", legend = c("CTL", "HFHC", "IRON", "IRON+HFHC"), col = c("#AFB9C8", "#8497B0","#F1CDB1","#F4B183"), ncol=2, 
       bty='n', cex = 0.8, pch = c(15, 16, 18, 17))
abline(v=c(20, 40), lty = 2)
dev.off()

################# FIGURE 1E ######################
colnames(meta_data[, c(2,3,5,11,12,18,23)])
res.pca <- PCA(meta_data[, c(2,3,5,18,23)], quali.sup = c(1,3), graph=T)

pdf(file = "figures/Figure_1E.pdf", width = 5, height = 5)
plot(res.pca$ind$coord, col = col.group, pch = pch.group, main = "PCA on variables describing steatosis",  cex.main=1)
abline(v=0, lty="dashed", col='grey')
abline(h=0, lty="dashed", col='grey')
legend("topright", legend = c("CTL", "HFHC", "IRON", "IRON+HFHC"), col = c("#AFB9C8", "#8497B0","#F1CDB1","#F4B183"), ncol=2, 
       bty='n', cex = 0.8, pch = c(15, 16, 18, 17))
dev.off()

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

pdf(file = "figures/Figure_1F.pdf", width = 8, height = 5)
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
              "Hepatic damages","Col1a mRNA","NRF2 mRNA","A2M mRNA","Crp mRNA","Total serum iron","Transferin Sat.","Hamp mRNA")
text(0.8:p,rep(-1.3, p), names.var, 
     srt = 45, xpd=NA, adj=c(1,1), cex = 0.8)
legend("topleft", legend = c("r correlation with Dim.1", "r correlation with Dim.2"),
       bty = 'n', pch=c(15, 16), cex=0.8, title.adj = 0)
text(x=c(1,6,9,11,14), y=-0.95, c("(1)", "(2)", "(3)", "(4)", "(5)"))
text(x=c(7.1), y=-0.95, "Used for\nthe PCA", cex=0.7)
dev.off()


