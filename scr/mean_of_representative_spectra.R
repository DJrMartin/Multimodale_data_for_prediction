rm(list=ls())
source("functions/preprocess.R")
#### representative spectra from the dataset #########
pheno <- read.csv("data/raw_data/phenotype.csv", sep=";", row.names=1)
df <- read.csv("data/raw_data/mir-frottis.csv", row.names = 1, sep=";")

# Time parameters.
TIME = df$Time
df = df[,-c(1)]
# There are several spectra from the same sample.
ind <- stringr::str_sub(rownames(df), end = -3)

# matplot(t(df), type="l")

df.mean = data.frame(matrix(NA, ncol = ncol(df), nrow=length(unique(ind))))
colnames(df.mean) = colnames(df)
for(i in 1:length(unique(ind))){
  if(sum(unique(ind)[i]==ind)>1){
    x <- df[unique(ind)[i]==ind,]
    # Computing the mean of the most representative spectra. 
    df.mean[i,] <- representative.mean(x, 0.2)$mean.spectra
  }else{
    # If there are not several spectra, we just get the one.
    df.mean[i,] <- df[i,]
  }
  rownames(df.mean)[1:i] = unique(ind)[1:i]
  print(i)
}

## Finalisation of the preprocessing
five.months <- df.mean[-grep("T5M", rownames(df.mean)),]
rownames(five.months) <- gsub("T8S","", rownames(five.months))
new.df <- five.months[rownames(pheno),]
rownames(new.df) = rownames(pheno)

# write.csv(new.df, file="data/data_clean/mir-frotti-mean-2months.csv")



