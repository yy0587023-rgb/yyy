###kinship
setwd("D:/gs2/10.16")
myY=read.table("AP.txt", head = TRUE)  
myY_hb  <- read.table("hb.txt",  head = TRUE, stringsAsFactors = FALSE)[,1:4]
myY_qy  <- read.table("qy.txt",  head = TRUE, stringsAsFactors = FALSE)[,1:4]
myY_hlj <- read.table("hlj.txt", head = TRUE, stringsAsFactors = FALSE)[,1:4]
myGD <- read.csv("myGD.csv", head = TRUE)
colnames(myGD)[1] <- "Taxa"
myY_hlj=myY_hlj[,1:4]
myY_hb=myY_hb[,1:4]
myY_hlj[1:10,]
myY_hb[1:10,]
remove_id<- c("D00592531", "D00608064")
myY_hb <- myY_hb[!myY_hb$Taxa %in% remove_id, ]
dim(myY_hb)
myY <- myY[!myY$Taxa %in% remove_id, ]
dim(myY)
remove_id1<- c("D00592537", "D00628064")
myGD<- myGD[!myGD$Taxa %in% remove_id1, ]
dim(myGD)
myGD_data <- myGD[, -1]
myGM <- read.csv("myGM.csv", head = TRUE)
myCV <- read.table("CV.txt", head = TRUE)[,1:3]
myCV<- myCV[!myCV$Taxa %in% remove_id, ]
dim(myCV)
source("D:/gs2/10.16/gapit_functions (6).txt") ##newest version
metal_path <- "D:/gs2/generic-metal/metal.exe"
nrep <- 1
nfold <- 277
cutOff <- 0.01
nSNP <- nrow(myGM)   
sig_thre <- cutOff / nSNP

gd_mat <- as.matrix(myGD[,-1])
gd_var <- apply(gd_mat, 2, var, na.rm=TRUE) 
gd_mat <- gd_mat[, gd_var > 0, drop=FALSE]
pca_res <- prcomp(gd_mat, scale. = TRUE)
PCs <- pca_res$x[,1:3]
colnames(PCs) <- paste0("PC", 1:3)
PCs_df <- data.frame(Taxa = myGD$Taxa, PCs)
head(PCs)

myCV1 <- data.frame(Taxa = myGD$Taxa, PCs)
head(myCV1)

PCs1<- data.frame(Taxa = myGD$Taxa, PCs)
myCV2<- merge(myCV, PCs1, by = "Taxa", all.x = TRUE)
head(myCV2)

results_all <- data.frame(
  trait = character(),
  rep = integer(),
  fold = integer(),
  method = character(),
  cor_pred = numeric(),
  cor_gBV  = numeric(),
  QTNs = integer(),
  stringsAsFactors = FALSE
)

traits <- c("DP")###
set.seed(1000)
for(trait in traits){
  cat("========== 正在预测性状:", trait, "==========\n")
    for(i in 1:nrep){

    sets_hb  <- sample(cut(1:nrow(myY_hb),  nfold, labels = FALSE), nrow(myY_hb))
    sets_qy  <- sample(cut(1:nrow(myY_qy),  nfold, labels = FALSE), nrow(myY_qy))
    sets_hlj <- sample(cut(1:nrow(myY_hlj), nfold, labels = FALSE), nrow(myY_hlj))
    pred_store_blink <- data.frame(Taxa=character(),
                                   Phenotype=numeric(),
                                   BLUE.N=numeric(),
                                   gBreedingValue=numeric(),
                                   stringsAsFactors = FALSE)

     pred_store_meta <- data.frame(Taxa=character(),
                                  Phenotype=numeric(),
                                  BLUE.N=numeric(),
                                  gBreedingValue=numeric(),
                                  stringsAsFactors = FALSE)
  
    for(j in 1:nfold){
      cat("rep=", i, " fold=", j, " trait=", trait, "\n")
     
      train_hb  <- myY_hb[ sets_hb  != j, c("Taxa", trait)]
      test_hb   <- myY_hb[ sets_hb  == j, c("Taxa", trait)]
      train_qy  <- myY_qy[ sets_qy  != j, c("Taxa", trait)]
      test_qy   <- myY_qy[ sets_qy  == j, c("Taxa", trait)]
      train_hlj <- myY_hlj[sets_hlj != j, c("Taxa", trait)]
      test_hlj  <- myY_hlj[sets_hlj == j, c("Taxa", trait)]
	  
      names(train_hb)[2]  <- "Phenotype"
      names(train_qy)[2]  <- "Phenotype"
      names(train_hlj)[2] <- "Phenotype"
      names(test_hb)[2]   <- "Phenotype"
      names(test_qy)[2]   <- "Phenotype"
      names(test_hlj)[2]  <- "Phenotype"
      
      # =================== (1)  GWAS ===================
      train_all <- rbind(train_hb, train_qy, train_hlj)
      #colnames(train_all) <- c("Taxa", "Phenotype")
      myGAPIT1 <- GAPIT(
        Y = train_all,
        GD = myGD,
        GM = myGM,
        CV = myCV2,  ##index+pca
        cutOff = cutOff,
        Random.model = FALSE,
        model = "Blink", 
		Multi_iter = FALSE
      )
      if (!is.null(myGAPIT1$GWAS) && "P.value" %in% names(myGAPIT1$GWAS)) {
      sig_res <- subset(myGAPIT1$GWAS, P.value < sig_thre)
      blink_sig_count <- nrow(sig_res)
      blink_sig_snps  <- intersect(sig_res$SNP, colnames(myGD))
}     else {
      blink_sig_count <- 0
      blink_sig_snps  <- character(0)
}

    cat("BLINK显著位点个数为:", blink_sig_count, "\n")
    if (length(blink_sig_snps) > 0) {
    myGD_sig_df <- data.frame(Taxa = myGD$Taxa, myGD[, blink_sig_snps, drop=FALSE])
    myCV3 <- merge(myCV1, myGD_sig_df, by = "Taxa", all.x = TRUE)
}   else {
    myCV3 <- myCV1
} 

      blink_gs <- GAPIT(
       Y = train_all,
       GD = myGD,
       GM = myGM,
       CV = myCV3, ##无index
       model=c("gBLUP")
)
      pred_blink <- blink_gs$Pred
	  test_all <- rbind(test_hb, test_qy, test_hlj)
      merged_blink <- merge(test_all, pred_blink[, c("Taxa", "BLUE.N", "gBreedingValue")],
                           by = "Taxa", all.x = TRUE)
	 pred_store_blink <- rbind(pred_store_blink,
                          merged_blink[, c("Taxa","Phenotype","BLUE.N","gBreedingValue")])

	  cor_pred_blink <- cor(merged_blink$Phenotype, merged_blink$BLUE.N)
      cor_gBV_blink<- cor(merged_blink$Phenotype, merged_blink$gBreedingValue)
      
	  cor_pred_blink
	  cor_gBV_blink
	  results_all <- rbind(results_all,
      data.frame(trait=trait, rep=i, fold=j, method="BLINK",
             cor_pred=cor_pred_blink, cor_gBV=cor_gBV_blink,
             QTNs = blink_sig_count,
             stringsAsFactors = FALSE))
      # =================== (2)meta+GS ===================
	  run_gwas <- function(train_df){
        res <- GAPIT(
          Y = train_df,
          GD = myGD,
          GM = myGM,
          CV = myCV1,
          cutOff = cutOff,
		  Random.model = FALSE,
          model = "Blink",
		  Multi_iter = FALSE
        )
        return(res$GWAS)
      }

      gwas_hb  <- run_gwas(train_hb)
      gwas_qy  <- run_gwas(train_qy)
      gwas_hlj <- run_gwas(train_hlj)

      # ======= 写入 METAL 输入文件 =======
      fn_hb  <- paste0("hb_rep", i, "_fold", j, ".txt")
      fn_qy  <- paste0("qy_rep", i, "_fold", j, ".txt")
      fn_hlj <- paste0("hlj_rep", i, "_fold", j, ".txt")

      write.table(gwas_hb[,c("SNP","effect","P.value","nobs")],  fn_hb,  sep="\t", row.names=FALSE, quote=FALSE)
      write.table(gwas_qy[,c("SNP","effect","P.value","nobs")],  fn_qy,  sep="\t", row.names=FALSE, quote=FALSE)
      write.table(gwas_hlj[,c("SNP","effect","P.value","nobs")], fn_hlj, sep="\t", row.names=FALSE, quote=FALSE)
      
      # ======= 写 METAL 脚本 =======
      metal_script <- paste0("metal_rep", i, "_fold", j, ".txt")
      meta_out <- paste0("meta_rep", i, "_fold", j)

      ms <- c(
        "SCHEME SAMPLESIZE",
        "MARKER SNP",
        "EFFECT effect",
        "PVALUE P.value",
        "WEIGHT nobs",
        paste0("PROCESS ", fn_hb),
        paste0("PROCESS ", fn_qy),
        paste0("PROCESS ", fn_hlj),
        paste0("OUTFILE ", meta_out),
        "ANALYZE",
        "QUIT"
      )
      writeLines(ms, metal_script)
    
	 system(paste(metal_path, metal_script), ignore.stdout=FALSE, ignore.stderr=FALSE)
     #taxa  <- myGD$taxa
     meta_tbl <- "METAANALYSIS1.TBL"
     meta_res <- read.table(meta_tbl, header=TRUE, sep="\t", stringsAsFactors=FALSE)
     sig_res  <- subset(meta_res, P.value < sig_thre)
     meta_sig_count <- nrow(sig_res)
     sig_snps <- intersect(sig_res$MarkerName, colnames(myGD_data))
 
     if (length(sig_snps) > 0) {
          myGD_sig_df <- data.frame(Taxa = myGD$Taxa, myGD_data[, sig_snps, drop = FALSE])
          myCV4 <- merge(myCV1, myGD_sig_df, by = "Taxa", all.x = TRUE)
      }
	else {
    myCV4 <- myCV1
} 

     meta_gs <- GAPIT(
       Y = train_all,
       GD = myGD,
       GM = myGM,
       CV = myCV4,
       model=c("gBLUP")
)
      pred_meta <- meta_gs$Pred
      merged_meta <- merge(test_all, pred_meta[, c("Taxa", "BLUE.N", "gBreedingValue")],
                           by = "Taxa", all.x = TRUE)

     pred_store_meta <- rbind(pred_store_meta,
                         merged_meta[, c("Taxa","Phenotype","BLUE.N","gBreedingValue")])
    if (j == 100) {
 
    out_file <- "D:/gs results/100kinship_gs.csv"
    if (!dir.exists(dirname(out_file))) dir.create(dirname(out_file), recursive = TRUE)
   
    pred_store_blink$method <- "BLINK"
    pred_store_meta$method  <- "META"
  
    pred_both <- rbind(
    pred_store_blink[, c("Taxa","Phenotype","BLUE.N","gBreedingValue","method")],
    pred_store_meta[,  c("Taxa","Phenotype","BLUE.N","gBreedingValue","method")]
  )
  
  write.csv(pred_both, out_file, row.names = FALSE)
  cat("✅ 已保存前 1..100 折的预测结果到文件：", out_file, "\n")
}
						   
      cor_pred_meta <- cor(merged_meta$Phenotype, merged_meta$BLUE.N)
      cor_gBV_meta  <- cor(merged_meta$Phenotype, merged_meta$gBreedingValue)
      
	  cor_pred_meta 
	  cor_gBV_meta
	  
      results_all <- rbind(results_all,
                           data.frame(trait=trait, rep=i, fold=j, method="Meta",
                                      cor_pred=cor_pred_meta, cor_gBV=cor_gBV_meta,
                                      QTNs = meta_sig_count, stringsAsFactors = FALSE))
									  
      
    cat(sprintf("Trait=%s | rep=%d | fold=%d | BLINK cor(Pred)=%.4f cor(gBV)=%.4f | Meta cor(Pred)=%.4f cor(gBV)=%.4f\n",
                  trait, i, j,
                  ifelse(is.na(cor_pred_blink), NaN, cor_pred_blink),
                  ifelse(is.na(cor_gBV_blink), NaN, cor_gBV_blink),
                  ifelse(is.na(cor_pred_meta), NaN, cor_pred_meta),
                  ifelse(is.na(cor_gBV_meta), NaN, cor_gBV_meta)))
}}}


if (!dir.exists("D:/gs results")) dir.create("D:/gs results", recursive = TRUE)
write.csv(results_all, "D:/gs results/6.0DP_blink_vs_meta(kinship).csv", row.names = FALSE)

##检查
setdiff(train_all$Taxa, myGD$Taxa)
setdiff(myY$Taxa, myGD$Taxa)
setdiff(myGD$Taxa, myY$Taxa) ##完全匹配


data<- read.csv("100kinship_gs.csv", head = TRUE)
data_blink <- subset(data, method == "BLINK")
dim(data_blink)
head(data_blink)
cor_pred_blink <- cor(data_blink$Phenotype, merged_blink$BLUE.N)
cor_gBV_blink <- cor(data_blink$Phenotype, merged_blink$gBreedingValue)

data_meta <- subset(data, method == "META")
dim(data_meta)
head(data_meta)
cor_pred_meta <- cor(data_blink$Phenotype, merged_meta$BLUE.N)
cor_gBV_meta <- cor(data_blink$Phenotype, merged_meta$gBreedingValue)