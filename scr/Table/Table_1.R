# renv::activate()
rm(list=ls())
library(FactoMineR)
library(dplyr)

load("data/data_preprocess/B-splines-normalisation.rda")
# load("data/data_preprocess/target_variables.rda")
# load("data/data_preprocess/fourier-normalisation.rda")
meta_data <- read.csv("data/raw_data/phenotype.csv", sep = ";")

metalomic <- read.csv("data/raw_data/metallomic.csv", sep = ",")
rownames(metalomic) = metalomic$X
metalomic <- metalomic[,-1]

df.explore <- meta_data[,c(2,9:17)] 
colnames(df$`metallomic-5months.csv`)[c(1,4,5,6,12:15)]
df.explore <- data.frame(Group = meta_data[,2], metalomic[meta_data$ID,c(1,4,5,6,12:15)])

# Chargement des packages nécessaires
library(dplyr)
library(tidyr)
library(ggpubr)
library(rstatix)

# Résumé statistique : Moyenne et écart-type par groupe
resume_stats <- df.explore %>%
  pivot_longer(cols = -Group, names_to = "Variable", values_to = "Valeur") %>%
  group_by(Group, Variable) %>%
  summarise(
    Moyenne = mean(Valeur, na.rm = TRUE),
    SD = sd(Valeur, na.rm = TRUE),
    .groups = "drop"
  )

# Fonction pour appliquer les tests statistiques par variable
dunn_results <- df.explore %>%
  pivot_longer(cols = -Group, names_to = "Variable", values_to = "Valeur") %>%
  group_by(Variable) %>%
  dunn_test(Valeur ~ Group, p.adjust.method = "bonferroni")

# Afficher les résultats
print("Résumé statistiques : Moyenne et Écart-type")
print(resume_stats)

print("Résultats des tests statistiques :")
print(resultats_stats)