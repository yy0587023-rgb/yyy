#GWAS
setwd("D:/gwas")
file_path ="gapit_functions.txt"
source(file_path)
install.packages("Biobase_2.66.0.zip", repos = NULL, type = "win.binary")
library(Biobase)

myY=read.table("AP.txt", head = TRUE)
myG=read.table("mdp.genotype.hmp.txt", head = FALSE)
myCV=read.csv("CV.csv", head = TRUE) 
myCV=myCV[,1:3]
myGAPIT <- GAPIT(
 G=myG,
 PCA.total=3,
 )
myGD= myGAPIT$GD
myGM= myGAPIT$GM
write.csv(myGD,"D:/gwas/myGD.csv",quote=F,row.names=F)
write.csv(myGM,"D:/gwas/myGM.csv",quote=F,row.names=F)

#shuffle-based method
myCV=read.csv("CV.csv", head = TRUE) 
myCV=myCV[,1:3]
Y <- myY[, c("Taxa", "CW")]  
n_perm <-50              
permuted_pvals <- numeric(n_perm)

for (i in 1:n_perm) {
  Y_permuted <- Y
  Y_permuted[, 2] <- sample(Y[, 2])  

  myGAPIT_perm <- GAPIT(
    Y = Y_permuted,
    GM = myGM,
    GD = myGD,
    CV = myCV,
    PCA.total = 3,
    model = c("BLINK"),
	)
  permuted_pvals[i] <- min(myGAPIT_perm$GWAS$P.value)
}

cat("All permutation min p-values:\n")
print(permuted_pvals)

alpha <- 0.05
threshold <- quantile(permuted_pvals, probs = alpha)
cat("Permutation-based shut-off threshold (50 times) for CW:", threshold, "\n")
#all populations
myGAPIT=GAPIT(
    Y =myY[,c(1,3)] ,
    GM = myGM,
    GD = myGD,
    CV = myCV,
    PCA.total = 3,
    model = c("BLINK"),
	cutOff=0.0861)
#single population
myY=read.table("hb.txt", head = TRUE)
myCV=read.csv("CV.csv", head = TRUE) 
myGD=read.csv("myGD.csv", head = TRUE) 
myGM=read.csv("myGM.csv", head = TRUE) 
myGAPIT=GAPIT(
  Y = myY[,c(1,3)],
  GM=myGM,
  GD=myGD,
  PCA.total = 3,
  model = c("BLINK"),
  cutOff=0.0861
 )
 
myY=read.table("qy.txt", head = TRUE)
myGAPIT=GAPIT(
  Y = myY[,c(1,5)],
  GM=myGM,
  GD=myGD, 
  PCA.total = 3,
  model = c("BLINK"),
  cutOff=0.0861
 )

myY=read.table("hlj.txt", head = TRUE)
myGAPIT=GAPIT(
  Y = myY[,c(1,3)],
  GM=myGM,
  GD=myGD,
  PCA.total = 3,
  model = c("BLINK"),
  cutOff=0.0861
 )

