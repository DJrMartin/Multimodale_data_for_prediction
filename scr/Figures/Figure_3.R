rm(list=ls())
library(yarrr)
library(caret)
library(randomForest)

load("data/results/LateIntegration.rda")
load("data/results/EarlyIntegration.rda")
load("data/results/MixedIntegration.rda")
load("data/results/SequentialIntegration.rda")

# pdf(file = "figures/Figure_3.pdf", width = 8.5, height = 4)
############### FIGURE 3A ###################
layout(matrix(c(1,1,2,3,4,5,5), nrow=1))
par(mar=c(10,4,3,1))
ylim = c(0.5, 2)

boxplot(data.frame(rmse.H0, rmse.blood, rmse.MIR, rmse.trace), axes=F, ylab = "RMSE", 
        col = c("grey", "#E1ECD8", "firebrick", "firebrick", "firebrick", ""), 
        ylim = ylim, main = "Models based on unimodal data", cex.main = 1)
box()
axis(1, at=1:4, labels=rep("",4))
text(1.3:4.3, 0.3, c("H0","Serum markers", "MIR (SelVar)", "Metallome"), 
     xpd = NA, cex = 0.8, srt = 45, pos = 2)
axis(2)
############### FIGURE 3B ###################
par(mar=c(10,1,3,1))
boxplot(data.frame(rmse.Early.XtremGB,rmse.Early.RF), 
        axes = F, ylab = "RMSE", col = c("white", "white"),
        ylim = ylim, main = "Early Integration", cex.main = 1)
box()
axis(1, at = 1:2, labels = rep("", 2))
text(1.3:3.2, 0.3, c("Xtrem GB","RF"), 
     xpd = NA, cex = 0.8, srt = 45, pos = 2)
############### FIGURE 3C ###################
par(mar=c(10,1,3,1))
boxplot(data.frame(rmse.Late), 
        axes = F, ylab = "RMSE", col = c("white"),
        ylim = ylim, main = "Late Integration", cex.main = 1)
box()
axis(1, at = 1:1, labels = rep("", 1))
text(1.2, 0.3, c("Linear integration"), 
     xpd = NA, cex = 0.8, srt = 45, pos = 2)
############### FIGURE 3D ###################
par(mar=c(10,1,3,1))
boxplot(data.frame(rmse.MixedR, rmse.MixedK), 
        axes = F, ylab = "RMSE", col = c("firebrick", "firebrick"),
        ylim = ylim, main = "Mixed Integration", cex.main = 1)
box()
axis(1, at = 1:2, labels = rep("", 2))
text(1.3:3.2, 0.3, c("Robust PCA","Kernel PCA"), 
     xpd = NA, cex = 0.8, srt = 45, pos = 2)
############### FIGURE 3E ###################
par(mar=c(10,1,3,1))
boxplot(data.frame(rmse.blood.MIR, rmse.blood.trace, rmse.trace.MIR, rmse.blood.MIR.trace, rmse.blood.trace.MIR), 
        axes=F, ylab = "RMSE", col = c("cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue"),
        ylim = ylim, main = "Sequential Integration", cex.main = 1)
t.test(rmse.blood, rmse.blood.trace.MIR)
box()
axis(1, at = 1:5, labels = rep("", 5))
text(1.3:5.3, 0.3, 
     c("Serum markers\n+ MIR (SelVar)", "Serum markers + metallome", "metallome + MIR (SelVar)", 
       "Serum markers + MIR (SelVar)\n+ metallome", "Serum markers + metallome\n+ MIR (SelVar)"),
     xpd = NA, cex = 0.8, srt = 45, pos = 2)

# dev.off()

# t.test(RMSE1, RMSE3)
pdf(file = "figures/legend_Fig_4.pdf", width = 8.5, height = 4.5)
plot.new()
legend("center", legend = c("Random" , "Neutral" ,"Decrease", "Increase"),title = "RMSE compared to RMSE from models fitted on Serum markers", 
       fill = c( "grey","white","cornflowerblue", "firebrick"), bty="n", ncol=4)
dev.off()