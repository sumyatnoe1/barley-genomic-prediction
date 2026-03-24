library(GAPIT)
library(readxl)

#### script_net blotch####

myY <- read.table("phenotypic_data.txt", head = T)
myG<- read.delim("genotypic_data.txt", header = F)


myGAPIT<- GAPIT(G=myG, output.numerical = TRUE)
myGD=myGAPIT$GD
myGM=myGAPIT$GM

#GAPIT run
myGAPIT_one <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  PCA.total = 3,SNP.MAF = 0.05,
  model=c("Blink"))