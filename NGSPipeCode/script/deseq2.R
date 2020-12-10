# DESeq2
# https://www.cnblogs.com/chenpeng1024/p/9260803.html

rm(list=ls())

library(DESeq2)
library(ggplot2)

# 读表达量矩阵
cts <- read.csv("/Users/zhangxuan/Work/Projects/bioinformatics/小桐子开花长非编码RNA/Analysis/RunAll/result/assembly_final/gene.tsv",sep="\t",row.names="Geneid")
colnames(cts) <- sub("\\.", "-", colnames(cts))

countData <- as.matrix(cts)
rownames(countData) <- rownames(cts)

# 读样本信息文件
coldata <- read.csv("/Users/zhangxuan/Work/Projects/bioinformatics/小桐子开花长非编码RNA/Analysis/RunAll/result/assembly_final/condition.tsv", row.names=1,sep="\t")
coldata <- coldata[,c("condition","type")]

all(rownames(coldata) %in% colnames(countData))
#countData <- countData[,rownames(coldata)]
all(rownames(coldata) == colnames(countData))

# 创建deseq2对象
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = coldata,
                              design = ~ condition)
dds

# 过滤低表达量的基因
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Differential expression analysis，差异表达分析的输入是原始reads文件，不能处理
## 运行deseq2
# 对原始dds进行标准化
dds <- DESeq(dds)
## 提取结果

diff_result_dir = "/Users/zhangxuan/Work/Projects/bioinformatics/小桐子开花长非编码RNA/Analysis/RunAll/result/diff123/"

sample_comb <- combn(unique(coldata$condition),2)

for (i in 1:dim(sample_comb)[2]){
  sample1 <- as.character(sample_comb[, i][1])
  sample2 <- as.character(sample_comb[, i][2])
  
  res <- results(dds, contrast=c("condition",sample1,sample2))
  
  resSig_all <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
  resSig_all_sorted <- resSig_all[ order(abs(resSig_all$log2FoldChange), decreasing = TRUE), ]
  write.csv(resSig_all_sorted,file= paste(diff_result_dir,sample1,"_vs_",sample2,"_all.csv", sep=''))
  
  resSig_up <- subset(res, padj < 0.05 & log2FoldChange > 1)
  resSig_up_sorted <- resSig_up[ order(resSig_up$log2FoldChange, decreasing = TRUE), ]
  write.csv(resSig_up_sorted,file= paste(diff_result_dir,sample1,"_vs_",sample2,"up.csv", sep=''))
  
  resSig_down <- subset(res, padj < 0.05 & log2FoldChange < 1)
  resSig_down_sorted <- resSig_down[ order(resSig_down$log2FoldChange), ]
  write.csv(resSig_down_sorted,file= paste(diff_result_dir,sample1,"_vs_",sample2,"down.csv", sep=''))
  
  # 火山图,在进行MA绘图之前，我们使用 lfcShrink函数缩小log2倍变化
  res.shrink <- lfcShrink(dds, contrast=c("condition",sample1,sample2), res=res)
  filename = paste(diff_result_dir,sample1,"_vs_",sample2,".pdf", sep='')
  cat(paste(sample1,"_vs_",sample2, sep=''))
  pdf(file=filename)
  plt <- plotMA(res.shrink, ylim = c(-5, 5))
  print(plt)
  dev.off()
}

# write.csv(resdata, "all_des_output.csv", row.names=FALSE)

