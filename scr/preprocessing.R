rm(list=ls())
source("functions/transform.R")

### Data normalisation with b-splines
path.to <- "data/data_clean"
df = list()
df[[dir(path.to)[1]]] = scale(read.csv( paste0(path.to ,"/", dir(path.to)[1]), row.names = 1))
for(i in dir(path.to)[-1]){
  df[[i]] = read.csv.transform(paste0(path.to,"/",i), m = 2, transformation = "splines", row.names = 1)
}
# save(df, file = "data/data_preprocess/b-splines-normalisation.rda")

