rm(list=ls())
library(yarrr)
library(caret)
library(randomForest)

load("data/results/LateIntegration.rda")
load("data/results/EarlyIntegration.rda")
load("data/results/MixedIntegration.rda")
load("data/results/SequentialIntegration.rda")

pdf(file = "figures/Figure_3.pdf", width = 8.5, height = 4)
############### FIGURE 3A ###################
layout(matrix(c(1,1,2,3,4,5,5), nrow=1))
par(mar=c(10,4,3,1))
ylim = c(0.5, 2)

boxplot(data.frame(rmse.H0, rmse.blood, rmse.MIR, rmse.trace), axes=F, ylab = "RMSE", 
        col = c("grey", "#E1ECD8", "firebrick", "white"), 
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
text(c(1,2), y=rep(1.8, 4), c("a,b", "a,b"))

############### FIGURE 3E ###################
par(mar=c(10,1,3,1))
boxplot(data.frame(rmse.blood.MIR, rmse.blood.trace, rmse.trace.MIR, rmse.blood.MIR.trace, rmse.blood.trace.MIR), 
        axes=F, ylab = "RMSE", col = c("cornflowerblue", "cornflowerblue", "white", "cornflowerblue", "cornflowerblue"),
        ylim = ylim, main = "Sequential Integration", cex.main = 1)
t.test(rmse.blood, rmse.blood.trace.MIR)
box()
axis(1, at = 1:5, labels = rep("", 5))
text(1.3:5.3, 0.3, 
     c("Serum markers\n+ MIR (SelVar)", "Serum markers + metallome", "metallome + MIR (SelVar)", 
       "Serum markers + MIR (SelVar)\n+ metallome", "Serum markers + metallome\n+ MIR (SelVar)"),
     xpd = NA, cex = 0.8, srt = 45, pos = 2)
text(c(1,2,4,5), y=rep(1.8, 4), c("b,c", "b,c", "b,c", "a,b,c"))
dev.off()

# pdf(file = "figures/legend_Fig_4.pdf", width = 8.5, height = 4.5)
# 
# plot.new()
# legend("center", legend = c("Random" , "Neutral" ,"Decrease", "Increase"),title = "RMSE compared to RMSE from models fitted on Serum markers", 
#        fill = c( "grey","white","cornflowerblue", "firebrick"), bty="n", ncol=4)
# 
# dev.off()

pairwise_ttest_bonferroni <- function(df, alpha = 0.05, paired = FALSE) {
  # Sélectionner uniquement les colonnes numériques
  numeric_df <- df[sapply(df, is.numeric)]
  col_names <- colnames(numeric_df)
  
  # Générer les combinaisons uniques de colonnes (sans doublons)
  combs <- combn(col_names, 2, simplify = FALSE)
  
  results <- lapply(combs, function(pair) {
    x <- numeric_df[[pair[1]]]
    y <- numeric_df[[pair[2]]]
    
    # Appliquer le t-test
    test <- try(t.test(x, y, paired = paired), silent = TRUE)
    
    if (!inherits(test, "try-error")) {
      return(data.frame(
        Var1 = pair[1],
        Var2 = pair[2],
        p_value = test$p.value
      ))
    } else {
      return(NULL)
    }
  })
  
  # Nettoyer et combiner les résultats
  results_df <- do.call(rbind, results)
  
  # Correction de Bonferroni
  results_df$p_adj <- p.adjust(results_df$p_value, method = "fdr")
  
  # Filtrer les tests significatifs
  signif_results <- subset(results_df, p_adj < alpha)
  
  return(signif_results)
}


df <- data.frame(rmse.blood, rmse.blood.MIR, rmse.blood.MIR.trace, rmse.blood.trace, rmse.blood.trace.MIR, 
                 rmse.Early.RF, rmse.Early.XtremGB, rmse.Late, rmse.MIR, rmse.MixedK, rmse.MixedR,
                 rmse.trace, rmse.trace.MIR)
pairwise_ttest_bonferroni(df, paired = T)


t.test(rmse.blood, rmse.blood.MIR.trace)

