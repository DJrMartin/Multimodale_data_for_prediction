renv::deactivate()
rm(list=ls())
library(FactoMineR)

# load("data/data_preprocess/B-splines-normalisation.rda")
load("data/data_preprocess/fourier-normalisation.rda")
meta_data <- read.csv("data/raw_data/phenotype.csv", sep = ";")

# matplot(t(df_2$`mir-frotti-mean-2months.csv`), type='l')

### STEATOSIS
colnames(meta_data[, c(2,3,11,12,18,23)])
res.pca <- PCA(meta_data[, c(2,3,18,23)], quali.sup = 1, graph=T)
# Dimension 1 is related to diet and steatosis score
# While Dim 2 is related to the iron supplementation.
synth.liver.steatosis <- res.pca$ind$coord[,1]
plot(meta_data$Triglycerides.Hp, meta_data$Steatosis)

### INFLAMMATION
colnames(meta_data[, c(2, 21, 22, 20, 24)])
res.pca <- PCA(meta_data[, c(2, 21, 22, 20, 24)], quali.sup = 1)
# Dim 1 is related to four transcris in liver which caracterize the liver inflammation.
synth.liver.inflammation <- res.pca$ind$coord[,1]
# Moreover this synthetic variable are significatively associated with the fibrosis in liver.
boxplot(synth.liver.inflammation~meta_data$Fibrose, xlab="Fibrosis", ylab=" Synthetic variable")
t.test(synth.liver.inflammation~meta_data$Fibrose)

### FINAL Synthetic variable 
res.pca <- PCA(data.frame(as.factor(meta_data$Group), as.factor(meta_data$Type), as.factor(meta_data$Fibrose), meta_data$Steatosis,
                     synth.liver.inflammation, synth.liver.steatosis), quali.sup = c(1:3), quanti.sup = 4)
# (1) We validate the synthetic variable for steatosis is well associated with histological score.
# Moreover this synthetic variable give an insight in the nature of the liver steatosis: 
# the mixte in macro and micro vesicules increase the severity of the steatosis. 
# (2) We also confirm that the synth.liver.inflammation is related to the fibrosis score.

S = data.frame(INFLAMMATION = synth.liver.inflammation, STEATOSIS = synth.liver.steatosis)
save(S ,file = "data/data_preprocess/target_variables.rda")

