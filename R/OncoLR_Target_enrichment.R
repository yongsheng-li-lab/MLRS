#' @title Enrichment analysis of LR inferred target genes based on clusterProfiler package.
#'
#' @param gene a vector of gene.
#' @param isRunGO Default TRUE. Run GO enrichment analysis.
#' @param isRunKEGG Default TRUE. Run KEGG enrichment analysis.
#'
#' @return
#' @export
#'
#' @examples
OncoLR_Target_enrichment <- function(gene,isRunGO=T,isRunKEGG=T){
  gene <- unique(gene)
  id_list <- mapIds(org.Hs.eg.db,gene,"ENTREZID","SYMBOL")
  id_list <- na.omit(id_list)
  if(isRunGO==T&isRunKEGG==F){
    print('Running GO enrichment')
    go <- enrichGO(gene = id_list, # Entrez ID列表
                   OrgDb = org.Hs.eg.db, # 指定物种数据库
                   keyType = "ENTREZID", # 指定给定的名称类型
                   ont = "ALL", # 可选,BP(生物学过程)/CC(细胞组分)/MF(分子功能)/ALL(同时指定)
                   pAdjustMethod = "BH", # P值校正方法,还可以是fdr
                   pvalueCutoff = 0.05,qvalueCutoff = 0.2, # p/q值阈值
                   readable = T # 将ID转换为symbol
    )
    go.res <- data.frame(go) # 将GO结果转为数据框，方便后续分析（不转也可以，看个人习惯）
    return(go.res)
  }
  if(isRunGO==F&isRunKEGG==T){
    print('Running KEGG enrichment')
    kegg <- enrichKEGG(gene = id_list,
                       organism = "hsa",keyType = "kegg",
                       pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2,
                       minGSSize = 10,maxGSSize = 500,use_internal_data = F)
    # 将结果表中的ID转换为symbol（非必需，根据个人需求）
    kk <- setReadable(kegg,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
    kegg.res <- data.frame(kk@result$ID,
                          kk@result$Description,
                          kk@result$GeneRatio,
                          kk@result$BgRatio,
                          kk@result$pvalue,
                          kk@result$p.adjust,
                          kk@result$geneID,
                          kk@result$Count)
    colnames(kegg.res) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","geneID","Count")
    return(kegg.res)
  }
  if(isRunGO==T&isRunKEGG==T){
    print('Step1----Running GO enrichment')
    go <- enrichGO(gene = id_list, # Entrez ID列表
                   OrgDb = org.Hs.eg.db, # 指定物种数据库
                   keyType = "ENTREZID", # 指定给定的名称类型
                   ont = "ALL", # 可选,BP(生物学过程)/CC(细胞组分)/MF(分子功能)/ALL(同时指定)
                   pAdjustMethod = "BH", # P值校正方法,还可以是fdr
                   pvalueCutoff = 0.05,qvalueCutoff = 0.2, # p/q值阈值
                   readable = T # 将ID转换为symbol
    )
    go.res <- data.frame(go) # 将GO结果转为数据框，方便后续分析（不转也可以，看个人习惯）
    print('Step2----Running KEGG enrichment')
    kegg <- enrichKEGG(gene = id_list,
                       organism = "hsa",keyType = "kegg",
                       pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2,
                       minGSSize = 10,maxGSSize = 500,use_internal_data = F)
    # 将结果表中的ID转换为symbol（非必需，根据个人需求）
    kk <- setReadable(kegg,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
    kegg.res <- data.frame(kk@result$ID,
                           kk@result$Description,
                           kk@result$GeneRatio,
                           kk@result$BgRatio,
                           kk@result$pvalue,
                           kk@result$p.adjust,
                           kk@result$geneID,
                           kk@result$Count)
    colnames(kegg.res) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","geneID","Count")
    return(list(go.res,kegg.res))
  }
}
