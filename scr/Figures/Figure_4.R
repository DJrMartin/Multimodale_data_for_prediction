rm(list=ls())
library(yarrr)
library(caret)
library(randomForest)

pdf(file = "figures/Figure_4.pdf", width = 8.5, height = 4)
############### FIGURE 4A ###################
layout(matrix(c(1,1,2,3,3), nrow=1))
par(mar=c(10,4,3,1))
ylim = c(range(res.A$rmse_per_fold, res.D.order$rmse_per_fold)[1], 2.2)

boxplot(data.frame(RMSE0, RMSE1, RMSE2, RMSE2.2, RMSE3, RMSE_all), axes=F, ylab = "RMSE", 
        col = c("grey", "white", "firebrick", "firebrick", "firebrick", "#E1ECD8"), 
        ylim = ylim, main = "Strategy n°1 - Pooled predictions", cex.main = 1)
box()
axis(1, at=1:6, labels=rep("",6))
text(1.3:6.3, 0.3, c("Null distribution","Serum markers","MIR","MIR (SelVar)","Metallome", "Pooled predictions"), 
     xpd = NA, cex = 0.8, srt = 45, pos = 2)
axis(2)
# segments(x0 = c(2,2,3), x1 = c(3,4,5), y0 = c(1.8, 1.7, 1.9))
# text(x = c(2.5, 3, 4), y = c(1.8, 1.7, 1.9) + 0.05, "*", cex = 1)

############### FIGURE 4B ###################
par(mar=c(10,1,3,2))
boxplot(data.frame(RMSE_all,  RMSE4, RMSE5), 
        axes = F, ylab = "RMSE", col = c("#E1ECD8", "white", "white"),
        ylim = ylim, main = "Strategy n°2 - Pooled data", cex.main = 1)
box()
abline(v = 1.5, lty = 2)
axis(1, at = 1:3, labels = rep("", 3))
text(1.3:3.3, 0.3, c("Serum markers","Gradient boosting","RF"), 
     xpd = NA, cex = 0.8, srt = 45, pos = 2)

############### FIGURE 4C ###################
par(mar=c(10,1,3,2))
boxplot(data.frame(RMSE_all,  res.A$rmse_per_fold, res.B$rmse_per_fold, 
                   res.C$rmse_per_fold, res.D$rmse_per_fold, res.D.order$rmse_per_fold), 
        axes=F, ylab = "RMSE", col = c("#E1ECD8", "cornflowerblue", "cornflowerblue", "white", "cornflowerblue", "cornflowerblue"),
        ylim = ylim, main = "Strategy n°3 - Aggregated models", cex.main = 1)
box()
abline(v = 1.5, lty = 'dashed')
axis(1, at = 1:6, labels = rep("", 6))
text(1.3:6.3, 0.3, 
     c("Serum markers","Serum markers\n+ MIR (SelVar)", "Serum markers + metallome", "MIR (SelVar) + metallome", 
       "Serum markers + MIR (SelVar)\n+ metallome", "Serum markers + metallome\n+ MIR (SelVar)"),
     xpd = NA, cex = 0.8, srt = 45, pos = 2)
# segments(x0 = c(5,1,1,3), x1 = c(6,5,3,6), y0= c(2,1.8, 1.7, 1.9))
# text(x = c(5.5), y=c(2.05), "***", cex=1)
# text(x = c(2), y=c(1.75), "*", cex=1)
# text(c(3, 4.5) , y = c(1.85, 1.95), "**", cex=1)
dev.off()

# t.test(RMSE1, RMSE3)
pdf(file = "figures/legend_Fig_4.pdf", width = 8.5, height = 4.5)
plot.new()
legend("center", legend = c("Random" , "Neutral" ,"Decrease", "Increase"),title = "RMSE compared to RMSE from models fitted on Serum markers", 
       fill = c( "grey","white","cornflowerblue", "firebrick"), bty="n", ncol=4)
dev.off()