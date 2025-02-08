rm(list = ls())

# prepare data
load('./lihc_tumor.Rdata')     # mRNA expression matrix (raw count) from TCGA-LIHC
load('./foxm1_target_trans_to_entrezID.Rdata')    # FOXM1 Target (gene symbols trans to entrezID)

# load data and packages
library(DESeq2)
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)

# build a results data frame
genes <- rownames(a)
gsea_res_df <- data.frame(genes)

# differential expression genes analysis and GSEA

for (i in 1:length(genes)){
  group_list <- ifelse(as.character(a[i,]>median(as.numeric(a[i,]))),
                       'high','low')
  group_list <- factor(group_list,levels = c('low','high'))
  
  colData <- data.frame(row.names = colnames(a), 
                        condition = group_list)
  dds <- DESeqDataSetFromMatrix(countData = a, 
                                colData = colData,
                                design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition",rev(levels(group_list))))
  resOrdered <- res[order(res$pvalue),]
  DEG <- as.data.frame(resOrdered)
  DEG <- na.omit(DEG)
  DEG <- DEG[!is.infinite(rowSums(DEG)),]
  
  gene <- bitr(rownames(DEG),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
  genes_df <- data.frame(SYMBOL= gene$SYMBOL,logFC=DEG[gene$SYMBOL,'log2FoldChange'])
  genes_df <- merge(genes_df,gene,by="SYMBOL")
  geneList <- genes_df$logFC
  names(geneList)=genes_df$ENTREZID
  geneList=sort(geneList,decreasing = T)
  

  tryCatch(
    {gsea_c <- GSEA(geneList = geneList, TERM2GENE = b,
                    verbose = FALSE,pvalueCutoff = 1) 
    gsea_res_df[i,2] <- gsea_c@result$NES
    gsea_res_df[i,3] <- gsea_c@result$pvalue},
    #   warning = function(w) {gsea_res_df[i,2:3] <- NA},
    error = function(e) {gsea_res_df[i,2:3] <- NA}
  )
  names(gsea_res_df) <- c("Symbol","NES","pvalue")
  if (i %% 1 == 0 ) {save(gsea_res_df,file = 'results.Rdata')
    
  } else {next
    
  }
  
}

# set the threshold and screen the genes 
pvalue_thres = 0.001
NES_thre = 1.5

gsea_res_sig <- gsea_res_df %>%
  filter(pvalue < pvalue_thres) %>% filter(abs(NES) > NES_thre) %>%
  arrange(desc(abs(NES)))
head(gsea_res_sig)
