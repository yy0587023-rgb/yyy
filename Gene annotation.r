#基因注释（100kb）
setwd("D:/gwas/基因注释")
snp_data=read.csv("significant markers for LW.csv",head=T)
#计算上下游100kb区间
snp_data$Start=snp_data$Pos - 50000
snp_data$End=snp_data$Pos + 50000
write.table(snp_data[, c("Chr", "Start", "End")], file="snp_regions_100.bed", sep="\t", row.names=FALSE, col.names=TRUE)
bed_data=read.table("snp_regions_100.bed", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(bed_data)

#读取GTF文件
#annotation=read.table("Bos_taurus.ARS-UCD2.0.114.gtf",
                         sep = "\t",
                         header = FALSE,
                         stringsAsFactors = FALSE,
                         fill = TRUE,
                         comment.char = "#")
#提取GTF中基因的信息
gene_annotation=annotation[annotation$V3 == "gene", ]
#从attribute字段中提取gene_id
gene_annotation$gene_id=sapply(gene_annotation$V9, function(x) {
  sub(".*gene_id \"([^\"]+)\".*", "\\1", x)
})
#提取感兴趣的列，构建成BED格式
#BED文件的start是从0开始计数，因此需要减去1
bed=data.frame(chrom = gene_annotation$V1,
                  start = gene_annotation$V4 - 1,
                  end = gene_annotation$V5,
                  gene_id = gene_annotation$gene_id)
write.table(bed, "gene_annotation.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#install.packages("data.table")
#install.packages("bedr")
library(data.table)
library(bedr)
snp_data1 <- fread("snp_regions_100.bed", header=TRUE)
snp_data1[, Start := as.integer(Start)]
snp_data1[, End := as.integer(End)]
head(snp_data1)
snp_data1[, Start := Start + 1]
gene_data1=fread("gene_annotation.bed", header=FALSE)
setnames(gene_data1, c("Chr", "Start", "End", "Info"))
head(gene_data1$Info)
#将 Info 列分解为单独的列
info_split=tstrsplit(gene_data1$Info, "; ", fixed=TRUE)
#创建新列并提取基因ID和基因名称
gene_data1[, ID := sub("gene_id ", "", info_split[[1]])]
gene_data1[, Name := sub("gene_name ", "", info_split[[3]])]
#删除 Info 列
gene_data1[, Info := NULL]
head(gene_data1)
gene_data1[, Start := Start + 1]
#定义函数来判断区间是否重叠
is_overlapping=function(snp_start, snp_end, gene_start, gene_end) {
  return(snp_start <= gene_end & snp_end >= gene_start)
}
#初始化结果数据框
result=list()
#遍历每个SNP区间，并检查是否与任何基因区间重叠
for (i in 1:nrow(snp_data1)) {
  snp=snp_data1[i, ]
  overlaps=gene_data1[
    is_overlapping(snp$Start, snp$End, Start, End),
  ]
  if (nrow(overlaps) > 0) {
    result[[length(result) + 1]]=cbind(snp, overlaps)
  }
}
library(data.table)
#合并结果
result_data=rbindlist(result)
#设置列名
setnames(result_data, c("SNP_Chr", "SNP_Start", "SNP_End",
                        "Gene_Chr", "Gene_Start", "Gene_End",
                        "Gene_ID", "Gene_Name"))
#过滤同染色体的数据
filtered_data=result_data[SNP_Chr == Gene_Chr]
#去重
unique_data=unique(filtered_data)
fwrite(unique_data, file = "LW_genes_50.bed", sep = "\t", quote = FALSE, col.names = TRUE)
