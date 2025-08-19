#' @title Calculate the effect of mutations on ligand/receptor functions.
#'
#' @param obj Defaults to ANNOVAR annotated results, the result file is a csv file.
#' @param LR_data A dataset with at least ligand/receptor columns of data. The default is the TICC and Omnipath datasets.
#' @param y_position1 numeric vector with the y positions of the brackets
#' @param ylim1 Limits for the y axes.
#' @param y_position2 numeric vector with the y positions of the brackets
#' @param ylim2 Limits for the y axes.
#' @param plot1 Default TRUE.
#' @param plot2 Default TRUE.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggpubr geom_signif
#'
#' @return
#' @export
#'
#' @examples
muteffect <- function(obj,LR_data,plot1=T,plot2=T,y_position1=NULL,ylim1=NULL,y_position2=NULL,ylim2=NULL){
  ##Reading data and organizing unique ligand/receptor
  ligand <- unique(L_R_data$ligand)
  receptor <- unique(L_R_data$receptor)
  ##Calculate Polyphen2 HDIV score
  cancer_mutation_Polyphen2_HDIV <- obj[,c("Hugo_Symbol","Polyphen2_HDIV_pred","Polyphen2_HDIV_score")]
  table(cancer_mutation_Polyphen2_HDIV$Polyphen2_HDIV_pred)
  cancer_mutation_Polyphen2_HDIV <- cancer_mutation_Polyphen2_HDIV[-which(cancer_mutation_Polyphen2_HDIV$Polyphen2_HDIV_pred=="."),]
  cancer_mutation_gene <- unique(cancer_mutation_Polyphen2_HDIV$Hugo_Symbol)
  cancer_mutation_ligand <- intersect(cancer_mutation_gene,ligand)
  cancer_mutation_ligand_other <- setdiff(cancer_mutation_gene,cancer_mutation_ligand)
  cancer_mutation_gene_L <- data.frame(c(cancer_mutation_ligand,cancer_mutation_ligand_other))
  colnames(cancer_mutation_gene_L) <- c("Hugo_Symbol")
  cancer_mutation_gene_L$L_lable <- c(rep("ligand",length(cancer_mutation_ligand)),
                                      rep("ligand_other",length(cancer_mutation_ligand_other)))
  cancer_mutation_Polyphen2_HDIV_L_LO <- merge(cancer_mutation_Polyphen2_HDIV,cancer_mutation_gene_L,by="Hugo_Symbol")
  cancer_mutation_receptor <- intersect(cancer_mutation_gene,receptor)
  cancer_mutation_receptor_other <- setdiff(cancer_mutation_gene,cancer_mutation_receptor)
  cancer_mutation_gene_R <- data.frame(c(cancer_mutation_receptor,cancer_mutation_receptor_other))
  colnames(cancer_mutation_gene_R) <- c("Hugo_Symbol")
  cancer_mutation_gene_R$R_lable <- c(rep("receptor",length(cancer_mutation_receptor)),
                                      rep("receptor_other",length(cancer_mutation_receptor_other)))
  cancer_mutation_Polyphen2_HDIV_R_RO <- merge(cancer_mutation_Polyphen2_HDIV,cancer_mutation_gene_R,by="Hugo_Symbol")
  colnames(cancer_mutation_Polyphen2_HDIV_L_LO) <- c(colnames(cancer_mutation_Polyphen2_HDIV_L_LO)[1:3],"lable")
  colnames(cancer_mutation_Polyphen2_HDIV_R_RO) <- c(colnames(cancer_mutation_Polyphen2_HDIV_R_RO)[1:3],"lable")
  cancer_Polyphen2_HDIV_box_data <- rbind(cancer_mutation_Polyphen2_HDIV_L_LO,
                                          cancer_mutation_Polyphen2_HDIV_R_RO)
  cancer_Polyphen2_HDIV_box_data$Polyphen2_HDIV_score <- as.numeric(cancer_Polyphen2_HDIV_box_data$Polyphen2_HDIV_score)
  colnames(cancer_Polyphen2_HDIV_box_data) <- c("Hugo_Symbol","Pred","Score",'lable')
  cancer_Polyphen2_HDIV_box_data$Method <- 'Polyphen2 HDIV'
  #wilcox.test
  cancer_mutation_Polyphen2_HDIV_L <- cancer_Polyphen2_HDIV_box_data[which(cancer_Polyphen2_HDIV_box_data$lable=='ligand'),]
  cancer_mutation_Polyphen2_HDIV_LO <- cancer_Polyphen2_HDIV_box_data[which(cancer_Polyphen2_HDIV_box_data$lable=='ligand_other'),]
  cancer_mutation_Polyphen2_HDIV_L_LO_P <- wilcox.test(cancer_mutation_Polyphen2_HDIV_L$Score,cancer_mutation_Polyphen2_HDIV_LO$Score)$p.value
  cancer_mutation_Polyphen2_HDIV_R <- cancer_Polyphen2_HDIV_box_data[which(cancer_Polyphen2_HDIV_box_data$lable=='receptor'),]
  cancer_mutation_Polyphen2_HDIV_RO <- cancer_Polyphen2_HDIV_box_data[which(cancer_Polyphen2_HDIV_box_data$lable=='receptor_other'),]
  cancer_mutation_Polyphen2_HDIV_R_RO_P <- wilcox.test(cancer_mutation_Polyphen2_HDIV_R$Score,cancer_mutation_Polyphen2_HDIV_RO$Score)$p.value
  ##Calculate Polyphen2 HVAR score
  cancer_mutation_Polyphen2_HVAR <- obj[,c("Hugo_Symbol","Polyphen2_HVAR_pred","Polyphen2_HVAR_score")]
  table(cancer_mutation_Polyphen2_HVAR$Polyphen2_HVAR_pred)
  cancer_mutation_Polyphen2_HVAR <- cancer_mutation_Polyphen2_HVAR[-which(cancer_mutation_Polyphen2_HVAR$Polyphen2_HVAR_pred=="."),]
  cancer_mutation_gene <- unique(cancer_mutation_Polyphen2_HVAR$Hugo_Symbol)
  cancer_mutation_ligand <- intersect(cancer_mutation_gene,ligand)
  cancer_mutation_ligand_other <- setdiff(cancer_mutation_gene,cancer_mutation_ligand)
  cancer_mutation_gene_L <- data.frame(c(cancer_mutation_ligand,cancer_mutation_ligand_other))
  colnames(cancer_mutation_gene_L) <- c("Hugo_Symbol")
  cancer_mutation_gene_L$L_lable <- c(rep("ligand",length(cancer_mutation_ligand)),
                                      rep("ligand_other",length(cancer_mutation_ligand_other)))
  cancer_mutation_Polyphen2_HVAR_L_LO <- merge(cancer_mutation_Polyphen2_HVAR,cancer_mutation_gene_L,by="Hugo_Symbol")
  cancer_mutation_receptor <- intersect(cancer_mutation_gene,receptor)
  cancer_mutation_receptor_other <- setdiff(cancer_mutation_gene,cancer_mutation_receptor)
  cancer_mutation_gene_R <- data.frame(c(cancer_mutation_receptor,cancer_mutation_receptor_other))
  colnames(cancer_mutation_gene_R) <- c("Hugo_Symbol")
  cancer_mutation_gene_R$R_lable <- c(rep("receptor",length(cancer_mutation_receptor)),
                                      rep("receptor_other",length(cancer_mutation_receptor_other)))
  cancer_mutation_Polyphen2_HVAR_R_RO <- merge(cancer_mutation_Polyphen2_HVAR,cancer_mutation_gene_R,by="Hugo_Symbol")
  colnames(cancer_mutation_Polyphen2_HVAR_L_LO) <- c(colnames(cancer_mutation_Polyphen2_HVAR_L_LO)[1:3],"lable")
  colnames(cancer_mutation_Polyphen2_HVAR_R_RO) <- c(colnames(cancer_mutation_Polyphen2_HVAR_R_RO)[1:3],"lable")
  cancer_Polyphen2_HVAR_box_data <- rbind(cancer_mutation_Polyphen2_HVAR_L_LO,
                                          cancer_mutation_Polyphen2_HVAR_R_RO)
  cancer_Polyphen2_HVAR_box_data$Polyphen2_HVAR_score <- as.numeric(cancer_Polyphen2_HVAR_box_data$Polyphen2_HVAR_score)
  colnames(cancer_Polyphen2_HVAR_box_data) <- c("Hugo_Symbol","Pred","Score",'lable')
  cancer_Polyphen2_HVAR_box_data$Method <- 'Polyphen2 HVAR'
  #wilcox.test
  cancer_mutation_Polyphen2_HVAR_L <- cancer_Polyphen2_HVAR_box_data[which(cancer_Polyphen2_HVAR_box_data$lable=='ligand'),]
  cancer_mutation_Polyphen2_HVAR_LO <- cancer_Polyphen2_HVAR_box_data[which(cancer_Polyphen2_HVAR_box_data$lable=='ligand_other'),]
  cancer_mutation_Polyphen2_HVAR_L_LO_P <- wilcox.test(cancer_mutation_Polyphen2_HVAR_L$Score,cancer_mutation_Polyphen2_HVAR_LO$Score)$p.value
  cancer_mutation_Polyphen2_HVAR_R <- cancer_Polyphen2_HVAR_box_data[which(cancer_Polyphen2_HVAR_box_data$lable=='receptor'),]
  cancer_mutation_Polyphen2_HVAR_RO <- cancer_Polyphen2_HVAR_box_data[which(cancer_Polyphen2_HVAR_box_data$lable=='receptor_other'),]
  cancer_mutation_Polyphen2_HVAR_R_RO_P <- wilcox.test(cancer_mutation_Polyphen2_HVAR_R$Score,cancer_mutation_Polyphen2_HVAR_RO$Score)$p.value
  ##Calculate LRT score
  cancer_mutation_LRT <- obj[,c("Hugo_Symbol","LRT_pred","LRT_score")]
  table(cancer_mutation_LRT$LRT_pred)
  cancer_mutation_LRT <- cancer_mutation_LRT[-which(cancer_mutation_LRT$LRT_pred=="."),]
  cancer_mutation_LRT <- cancer_mutation_LRT[-which(cancer_mutation_LRT$LRT_pred=="U"),]
  cancer_mutation_gene <- unique(cancer_mutation_LRT$Hugo_Symbol)
  cancer_mutation_ligand <- intersect(cancer_mutation_gene,ligand)
  cancer_mutation_ligand_other <- setdiff(cancer_mutation_gene,cancer_mutation_ligand)
  cancer_mutation_gene_L <- data.frame(c(cancer_mutation_ligand,cancer_mutation_ligand_other))
  colnames(cancer_mutation_gene_L) <- c("Hugo_Symbol")
  cancer_mutation_gene_L$L_lable <- c(rep("ligand",length(cancer_mutation_ligand)),
                                      rep("ligand_other",length(cancer_mutation_ligand_other)))
  cancer_mutation_LRT_L_LO <- merge(cancer_mutation_LRT,cancer_mutation_gene_L,by="Hugo_Symbol")
  cancer_mutation_receptor <- intersect(cancer_mutation_gene,receptor)
  cancer_mutation_receptor_other <- setdiff(cancer_mutation_gene,cancer_mutation_receptor)
  cancer_mutation_gene_R <- data.frame(c(cancer_mutation_receptor,cancer_mutation_receptor_other))
  colnames(cancer_mutation_gene_R) <- c("Hugo_Symbol")
  cancer_mutation_gene_R$R_lable <- c(rep("receptor",length(cancer_mutation_receptor)),
                                      rep("receptor_other",length(cancer_mutation_receptor_other)))
  cancer_mutation_LRT_R_RO <- merge(cancer_mutation_LRT,cancer_mutation_gene_R,by="Hugo_Symbol")
  colnames(cancer_mutation_LRT_L_LO) <- c(colnames(cancer_mutation_LRT_L_LO)[1:3],"lable")
  colnames(cancer_mutation_LRT_R_RO) <- c(colnames(cancer_mutation_LRT_R_RO)[1:3],"lable")
  cancer_LRT_box_data <- rbind(cancer_mutation_LRT_L_LO,
                               cancer_mutation_LRT_R_RO)
  cancer_LRT_box_data$LRT_score <- as.numeric(cancer_LRT_box_data$LRT_score)
  colnames(cancer_LRT_box_data) <- c("Hugo_Symbol","Pred","Score",'lable')
  cancer_LRT_box_data$Method <- 'LRT'
  #wilcox.test
  cancer_mutation_LRT_L <- cancer_LRT_box_data[which(cancer_LRT_box_data$lable=='ligand'),]
  cancer_mutation_LRT_LO <- cancer_LRT_box_data[which(cancer_LRT_box_data$lable=='ligand_other'),]
  cancer_mutation_LRT_L_LO_P <- wilcox.test(cancer_mutation_LRT_L$Score,cancer_mutation_LRT_LO$Score)$p.value
  cancer_mutation_LRT_R <- cancer_LRT_box_data[which(cancer_LRT_box_data$lable=='receptor'),]
  cancer_mutation_LRT_RO <- cancer_LRT_box_data[which(cancer_LRT_box_data$lable=='receptor_other'),]
  cancer_mutation_LRT_R_RO_P <- wilcox.test(cancer_mutation_LRT_R$Score,cancer_mutation_LRT_RO$Score)$p.value
  ##Calculate MetaSVM score
  cancer_mutation_MetaSVM <- obj[,c("Hugo_Symbol","MetaSVM_pred","MetaSVM_score")]
  table(cancer_mutation_MetaSVM$MetaSVM_pred)
  cancer_mutation_MetaSVM <- cancer_mutation_MetaSVM[-which(cancer_mutation_MetaSVM$MetaSVM_pred=="."),]
  cancer_mutation_gene <- unique(cancer_mutation_MetaSVM$Hugo_Symbol)
  cancer_mutation_ligand <- intersect(cancer_mutation_gene,ligand)
  cancer_mutation_ligand_other <- setdiff(cancer_mutation_gene,cancer_mutation_ligand)
  cancer_mutation_gene_L <- data.frame(c(cancer_mutation_ligand,cancer_mutation_ligand_other))
  colnames(cancer_mutation_gene_L) <- c("Hugo_Symbol")
  cancer_mutation_gene_L$L_lable <- c(rep("ligand",length(cancer_mutation_ligand)),
                                      rep("ligand_other",length(cancer_mutation_ligand_other)))
  cancer_mutation_MetaSVM_L_LO <- merge(cancer_mutation_MetaSVM,cancer_mutation_gene_L,by="Hugo_Symbol")
  cancer_mutation_receptor <- intersect(cancer_mutation_gene,receptor)
  cancer_mutation_receptor_other <- setdiff(cancer_mutation_gene,cancer_mutation_receptor)
  cancer_mutation_gene_R <- data.frame(c(cancer_mutation_receptor,cancer_mutation_receptor_other))
  colnames(cancer_mutation_gene_R) <- c("Hugo_Symbol")
  cancer_mutation_gene_R$R_lable <- c(rep("receptor",length(cancer_mutation_receptor)),
                                      rep("receptor_other",length(cancer_mutation_receptor_other)))
  cancer_mutation_MetaSVM_R_RO <- merge(cancer_mutation_MetaSVM,cancer_mutation_gene_R,by="Hugo_Symbol")
  colnames(cancer_mutation_MetaSVM_L_LO) <- c(colnames(cancer_mutation_MetaSVM_L_LO)[1:3],"lable")
  colnames(cancer_mutation_MetaSVM_R_RO) <- c(colnames(cancer_mutation_MetaSVM_R_RO)[1:3],"lable")
  cancer_MetaSVM_box_data <- rbind(cancer_mutation_MetaSVM_L_LO,
                                   cancer_mutation_MetaSVM_R_RO)
  cancer_MetaSVM_box_data$MetaSVM_score <- as.numeric(cancer_MetaSVM_box_data$MetaSVM_score)
  colnames(cancer_MetaSVM_box_data) <- c("Hugo_Symbol","Pred","Score",'lable')
  cancer_MetaSVM_box_data$Method <- 'MetaSVM'
  #wilcox.test
  cancer_mutation_MetaSVM_L <- cancer_MetaSVM_box_data[which(cancer_MetaSVM_box_data$lable=='ligand'),]
  cancer_mutation_MetaSVM_LO <- cancer_MetaSVM_box_data[which(cancer_MetaSVM_box_data$lable=='ligand_other'),]
  cancer_mutation_MetaSVM_L_LO_P <- wilcox.test(cancer_mutation_MetaSVM_L$Score,cancer_mutation_MetaSVM_LO$Score)$p.value
  cancer_mutation_MetaSVM_R <- cancer_MetaSVM_box_data[which(cancer_MetaSVM_box_data$lable=='receptor'),]
  cancer_mutation_MetaSVM_RO <- cancer_MetaSVM_box_data[which(cancer_MetaSVM_box_data$lable=='receptor_other'),]
  cancer_mutation_MetaSVM_R_RO_P <- wilcox.test(cancer_mutation_MetaSVM_R$Score,cancer_mutation_MetaSVM_RO$Score)$p.value
  ##Calculate MetaLR score
  cancer_mutation_MetaLR <- obj[,c("Hugo_Symbol","MetaLR_pred","MetaLR_score")]
  table(cancer_mutation_MetaLR$MetaLR_pred)
  cancer_mutation_MetaLR <- cancer_mutation_MetaLR[-which(cancer_mutation_MetaLR$MetaLR_pred=="."),]
  cancer_mutation_gene <- unique(cancer_mutation_MetaLR$Hugo_Symbol)
  cancer_mutation_ligand <- intersect(cancer_mutation_gene,ligand)
  cancer_mutation_ligand_other <- setdiff(cancer_mutation_gene,cancer_mutation_ligand)
  cancer_mutation_gene_L <- data.frame(c(cancer_mutation_ligand,cancer_mutation_ligand_other))
  colnames(cancer_mutation_gene_L) <- c("Hugo_Symbol")
  cancer_mutation_gene_L$L_lable <- c(rep("ligand",length(cancer_mutation_ligand)),
                                      rep("ligand_other",length(cancer_mutation_ligand_other)))
  cancer_mutation_MetaLR_L_LO <- merge(cancer_mutation_MetaLR,cancer_mutation_gene_L,by="Hugo_Symbol")
  cancer_mutation_receptor <- intersect(cancer_mutation_gene,receptor)
  cancer_mutation_receptor_other <- setdiff(cancer_mutation_gene,cancer_mutation_receptor)
  cancer_mutation_gene_R <- data.frame(c(cancer_mutation_receptor,cancer_mutation_receptor_other))
  colnames(cancer_mutation_gene_R) <- c("Hugo_Symbol")
  cancer_mutation_gene_R$R_lable <- c(rep("receptor",length(cancer_mutation_receptor)),
                                      rep("receptor_other",length(cancer_mutation_receptor_other)))
  cancer_mutation_MetaLR_R_RO <- merge(cancer_mutation_MetaLR,cancer_mutation_gene_R,by="Hugo_Symbol")
  colnames(cancer_mutation_MetaLR_L_LO) <- c(colnames(cancer_mutation_MetaLR_L_LO)[1:3],"lable")
  colnames(cancer_mutation_MetaLR_R_RO) <- c(colnames(cancer_mutation_MetaLR_R_RO)[1:3],"lable")
  cancer_MetaLR_box_data <- rbind(cancer_mutation_MetaLR_L_LO,
                                  cancer_mutation_MetaLR_R_RO)
  cancer_MetaLR_box_data$MetaLR_score <- as.numeric(cancer_MetaLR_box_data$MetaLR_score)
  colnames(cancer_MetaLR_box_data) <- c("Hugo_Symbol","Pred","Score",'lable')
  cancer_MetaLR_box_data$Method <- 'MetaLR'
  #wilcox.test
  cancer_mutation_MetaLR_L <- cancer_MetaLR_box_data[which(cancer_MetaLR_box_data$lable=='ligand'),]
  cancer_mutation_MetaLR_LO <- cancer_MetaLR_box_data[which(cancer_MetaLR_box_data$lable=='ligand_other'),]
  cancer_mutation_MetaLR_L_LO_P <- wilcox.test(cancer_mutation_MetaLR_L$Score,cancer_mutation_MetaLR_LO$Score)$p.value
  cancer_mutation_MetaLR_R <- cancer_MetaLR_box_data[which(cancer_MetaLR_box_data$lable=='receptor'),]
  cancer_mutation_MetaLR_RO <- cancer_MetaLR_box_data[which(cancer_MetaLR_box_data$lable=='receptor_other'),]
  cancer_mutation_MetaLR_R_RO_P <- wilcox.test(cancer_mutation_MetaLR_R$Score,cancer_mutation_MetaLR_RO$Score)$p.value
  ##Calculate M CAP score
  cancer_mutation_M_CAP <- obj[,c("Hugo_Symbol","M.CAP_pred","M.CAP_score")]
  colnames(cancer_mutation_M_CAP) <- c(c("Hugo_Symbol","M_CAP_pred","M_CAP_score"))
  table(cancer_mutation_M_CAP$M_CAP_pred)
  cancer_mutation_M_CAP <- cancer_mutation_M_CAP[-which(cancer_mutation_M_CAP$M_CAP_pred=="."),]
  cancer_mutation_gene <- unique(cancer_mutation_M_CAP$Hugo_Symbol)
  cancer_mutation_ligand <- intersect(cancer_mutation_gene,ligand)
  cancer_mutation_ligand_other <- setdiff(cancer_mutation_gene,cancer_mutation_ligand)
  cancer_mutation_gene_L <- data.frame(c(cancer_mutation_ligand,cancer_mutation_ligand_other))
  colnames(cancer_mutation_gene_L) <- c("Hugo_Symbol")
  cancer_mutation_gene_L$L_lable <- c(rep("ligand",length(cancer_mutation_ligand)),
                                      rep("ligand_other",length(cancer_mutation_ligand_other)))
  cancer_mutation_M_CAP_L_LO <- merge(cancer_mutation_M_CAP,cancer_mutation_gene_L,by="Hugo_Symbol")
  cancer_mutation_receptor <- intersect(cancer_mutation_gene,receptor)
  cancer_mutation_receptor_other <- setdiff(cancer_mutation_gene,cancer_mutation_receptor)
  cancer_mutation_gene_R <- data.frame(c(cancer_mutation_receptor,cancer_mutation_receptor_other))
  colnames(cancer_mutation_gene_R) <- c("Hugo_Symbol")
  cancer_mutation_gene_R$R_lable <- c(rep("receptor",length(cancer_mutation_receptor)),
                                      rep("receptor_other",length(cancer_mutation_receptor_other)))
  cancer_mutation_M_CAP_R_RO <- merge(cancer_mutation_M_CAP,cancer_mutation_gene_R,by="Hugo_Symbol")
  colnames(cancer_mutation_M_CAP_L_LO) <- c(colnames(cancer_mutation_M_CAP_L_LO)[1:3],"lable")
  colnames(cancer_mutation_M_CAP_R_RO) <- c(colnames(cancer_mutation_M_CAP_R_RO)[1:3],"lable")
  cancer_M_CAP_box_data <- rbind(cancer_mutation_M_CAP_L_LO,
                                 cancer_mutation_M_CAP_R_RO)
  cancer_M_CAP_box_data$M_CAP_score <- as.numeric(cancer_M_CAP_box_data$M_CAP_score)
  colnames(cancer_M_CAP_box_data) <- c("Hugo_Symbol","Pred","Score",'lable')
  cancer_M_CAP_box_data$Method <- 'M CAP'
  #wilcox.test
  cancer_mutation_M_CAP_L <- cancer_M_CAP_box_data[which(cancer_M_CAP_box_data$lable=='ligand'),]
  cancer_mutation_M_CAP_LO <- cancer_M_CAP_box_data[which(cancer_M_CAP_box_data$lable=='ligand_other'),]
  cancer_mutation_M_CAP_L_LO_P <- wilcox.test(cancer_mutation_M_CAP_L$Score,cancer_mutation_M_CAP_LO$Score)$p.value
  cancer_mutation_M_CAP_R <- cancer_M_CAP_box_data[which(cancer_M_CAP_box_data$lable=='receptor'),]
  cancer_mutation_M_CAP_RO <- cancer_M_CAP_box_data[which(cancer_M_CAP_box_data$lable=='receptor_other'),]
  cancer_mutation_M_CAP_R_RO_P <- wilcox.test(cancer_mutation_M_CAP_R$Score,cancer_mutation_M_CAP_RO$Score)$p.value
  ##Calculate DEOGEN2 score
  cancer_mutation_DEOGEN2 <- obj[,c("Hugo_Symbol","DEOGEN2_pred","DEOGEN2_score")]
  table(cancer_mutation_DEOGEN2$DEOGEN2_pred)
  cancer_mutation_DEOGEN2 <- cancer_mutation_DEOGEN2[-which(cancer_mutation_DEOGEN2$DEOGEN2_pred=="."),]
  cancer_mutation_gene <- unique(cancer_mutation_DEOGEN2$Hugo_Symbol)
  cancer_mutation_ligand <- intersect(cancer_mutation_gene,ligand)
  cancer_mutation_ligand_other <- setdiff(cancer_mutation_gene,cancer_mutation_ligand)
  cancer_mutation_gene_L <- data.frame(c(cancer_mutation_ligand,cancer_mutation_ligand_other))
  colnames(cancer_mutation_gene_L) <- c("Hugo_Symbol")
  cancer_mutation_gene_L$L_lable <- c(rep("ligand",length(cancer_mutation_ligand)),
                                      rep("ligand_other",length(cancer_mutation_ligand_other)))
  cancer_mutation_DEOGEN2_L_LO <- merge(cancer_mutation_DEOGEN2,cancer_mutation_gene_L,by="Hugo_Symbol")
  cancer_mutation_receptor <- intersect(cancer_mutation_gene,receptor)
  cancer_mutation_receptor_other <- setdiff(cancer_mutation_gene,cancer_mutation_receptor)
  cancer_mutation_gene_R <- data.frame(c(cancer_mutation_receptor,cancer_mutation_receptor_other))
  colnames(cancer_mutation_gene_R) <- c("Hugo_Symbol")
  cancer_mutation_gene_R$R_lable <- c(rep("receptor",length(cancer_mutation_receptor)),
                                      rep("receptor_other",length(cancer_mutation_receptor_other)))
  cancer_mutation_DEOGEN2_R_RO <- merge(cancer_mutation_DEOGEN2,cancer_mutation_gene_R,by="Hugo_Symbol")
  colnames(cancer_mutation_DEOGEN2_L_LO) <- c(colnames(cancer_mutation_DEOGEN2_L_LO)[1:3],"lable")
  colnames(cancer_mutation_DEOGEN2_R_RO) <- c(colnames(cancer_mutation_DEOGEN2_R_RO)[1:3],"lable")
  cancer_DEOGEN2_box_data <- rbind(cancer_mutation_DEOGEN2_L_LO,
                                   cancer_mutation_DEOGEN2_R_RO)
  cancer_DEOGEN2_box_data$DEOGEN2_score <- as.numeric(cancer_DEOGEN2_box_data$DEOGEN2_score)
  colnames(cancer_DEOGEN2_box_data) <- c("Hugo_Symbol","Pred","Score",'lable')
  cancer_DEOGEN2_box_data$Method <- 'DEOGEN2'
  #wilcox.test
  cancer_mutation_DEOGEN2_L <- cancer_DEOGEN2_box_data[which(cancer_DEOGEN2_box_data$lable=='ligand'),]
  cancer_mutation_DEOGEN2_LO <- cancer_DEOGEN2_box_data[which(cancer_DEOGEN2_box_data$lable=='ligand_other'),]
  cancer_mutation_DEOGEN2_L_LO_P <- wilcox.test(cancer_mutation_DEOGEN2_L$Score,cancer_mutation_DEOGEN2_LO$Score)$p.value
  cancer_mutation_DEOGEN2_R <- cancer_DEOGEN2_box_data[which(cancer_DEOGEN2_box_data$lable=='receptor'),]
  cancer_mutation_DEOGEN2_RO <- cancer_DEOGEN2_box_data[which(cancer_DEOGEN2_box_data$lable=='receptor_other'),]
  cancer_mutation_DEOGEN2_R_RO_P <- wilcox.test(cancer_mutation_DEOGEN2_R$Score,cancer_mutation_DEOGEN2_RO$Score)$p.value
  #Integrate damage_data
  cancer_damage_data <- rbind(cancer_Polyphen2_HDIV_box_data,
                              cancer_Polyphen2_HVAR_box_data,
                              cancer_LRT_box_data,
                              cancer_MetaLR_box_data,
                              cancer_M_CAP_box_data,
                              cancer_DEOGEN2_box_data)
  cancer_damage_data <- cancer_damage_data[,-2]
  cancer_mutation_method_wilcox_P <- data.frame(c(cancer_mutation_Polyphen2_HDIV_L_LO_P,
                                                  cancer_mutation_Polyphen2_HDIV_R_RO_P,
                                                  cancer_mutation_Polyphen2_HVAR_L_LO_P,
                                                  cancer_mutation_Polyphen2_HVAR_R_RO_P,
                                                  cancer_mutation_DEOGEN2_L_LO_P,
                                                  cancer_mutation_DEOGEN2_R_RO_P,
                                                  cancer_mutation_LRT_L_LO_P,
                                                  cancer_mutation_LRT_R_RO_P,
                                                  cancer_mutation_M_CAP_L_LO_P,
                                                  cancer_mutation_M_CAP_R_RO_P,
                                                  cancer_mutation_MetaLR_L_LO_P,
                                                  cancer_mutation_MetaLR_R_RO_P
  ))
  colnames(cancer_mutation_method_wilcox_P) <- c("P")
  cancer_mutation_method_wilcox_P$Method <- c('Polyphen2_HDIV_L_LO',
                                              'Polyphen2_HDIV_R_RO',
                                              'Polyphen2_HVAR_L_LO',
                                              'Polyphen2_HVAR_R_RO',
                                              'DEOGEN2_L_LO',
                                              'DEOGEN2_R_RO',
                                              'LRT_L_LO',
                                              'LRT_R_RO',
                                              'M_CAP_L_LO',
                                              'M_CAP_R_RO',
                                              'MetaLR_L_LO',
                                              'MetaLR_R_RO')
  cancer_mutation_method_wilcox_P$p_lable <- ifelse(cancer_mutation_method_wilcox_P$P>0.05,"NS",
                                                    ifelse(cancer_mutation_method_wilcox_P$P>0.01,"*",
                                                           ifelse(cancer_mutation_method_wilcox_P$P>0.001,"**",
                                                                  ifelse(cancer_mutation_method_wilcox_P$P>0.0001,"***","****"))))
  ##Calculate phastCons100way vertebrate score
  cancer_mutation_phastCons100way_vertebrate <- obj[,c("Hugo_Symbol","phastCons100way_vertebrate")]
  colnames(cancer_mutation_phastCons100way_vertebrate) <- c("Hugo_Symbol","phastCons100way_vertebrate_score")
  cancer_mutation_phastCons100way_vertebrate <- cancer_mutation_phastCons100way_vertebrate[-which(cancer_mutation_phastCons100way_vertebrate$phastCons100way_vertebrate_score=="."),]
  cancer_mutation_gene <- unique(cancer_mutation_phastCons100way_vertebrate$Hugo_Symbol)
  cancer_mutation_ligand <- intersect(cancer_mutation_gene,ligand)
  cancer_mutation_ligand_other <- setdiff(cancer_mutation_gene,cancer_mutation_ligand)
  cancer_mutation_gene_L <- data.frame(c(cancer_mutation_ligand,cancer_mutation_ligand_other))
  colnames(cancer_mutation_gene_L) <- c("Hugo_Symbol")
  cancer_mutation_gene_L$L_lable <- c(rep("ligand",length(cancer_mutation_ligand)),
                                      rep("ligand_other",length(cancer_mutation_ligand_other)))
  cancer_mutation_phastCons100way_vertebrate_L_LO <- merge(cancer_mutation_phastCons100way_vertebrate,cancer_mutation_gene_L,by="Hugo_Symbol")
  cancer_mutation_receptor <- intersect(cancer_mutation_gene,receptor)
  cancer_mutation_receptor_other <- setdiff(cancer_mutation_gene,cancer_mutation_receptor)
  cancer_mutation_gene_R <- data.frame(c(cancer_mutation_receptor,cancer_mutation_receptor_other))
  colnames(cancer_mutation_gene_R) <- c("Hugo_Symbol")
  cancer_mutation_gene_R$R_lable <- c(rep("receptor",length(cancer_mutation_receptor)),
                                      rep("receptor_other",length(cancer_mutation_receptor_other)))
  cancer_mutation_phastCons100way_vertebrate_R_RO <- merge(cancer_mutation_phastCons100way_vertebrate,cancer_mutation_gene_R,by="Hugo_Symbol")
  colnames(cancer_mutation_phastCons100way_vertebrate_L_LO) <- c(colnames(cancer_mutation_phastCons100way_vertebrate_L_LO)[1:2],"lable")
  colnames(cancer_mutation_phastCons100way_vertebrate_R_RO) <- c(colnames(cancer_mutation_phastCons100way_vertebrate_R_RO)[1:2],"lable")
  cancer_phastCons100way_vertebrate_box_data <- rbind(cancer_mutation_phastCons100way_vertebrate_L_LO,
                                                      cancer_mutation_phastCons100way_vertebrate_R_RO)
  cancer_phastCons100way_vertebrate_box_data$phastCons100way_vertebrate_score <- as.numeric(cancer_phastCons100way_vertebrate_box_data$phastCons100way_vertebrate_score)
  colnames(cancer_phastCons100way_vertebrate_box_data) <- c("Hugo_Symbol","Score",'lable')
  cancer_phastCons100way_vertebrate_box_data$Method <- 'phastCons100way vertebrate'
  #wilcox.test
  cancer_mutation_phastCons100way_vertebrate_L <- cancer_phastCons100way_vertebrate_box_data[which(cancer_phastCons100way_vertebrate_box_data$lable=='ligand'),]
  cancer_mutation_phastCons100way_vertebrate_LO <- cancer_phastCons100way_vertebrate_box_data[which(cancer_phastCons100way_vertebrate_box_data$lable=='ligand_other'),]
  cancer_mutation_phastCons100way_vertebrate_L_LO_P <- wilcox.test(cancer_mutation_phastCons100way_vertebrate_L$Score,cancer_mutation_phastCons100way_vertebrate_LO$Score)$p.value
  cancer_mutation_phastCons100way_vertebrate_R <- cancer_phastCons100way_vertebrate_box_data[which(cancer_phastCons100way_vertebrate_box_data$lable=='receptor'),]
  cancer_mutation_phastCons100way_vertebrate_RO <- cancer_phastCons100way_vertebrate_box_data[which(cancer_phastCons100way_vertebrate_box_data$lable=='receptor_other'),]
  cancer_mutation_phastCons100way_vertebrate_R_RO_P <- wilcox.test(cancer_mutation_phastCons100way_vertebrate_R$Score,cancer_mutation_phastCons100way_vertebrate_RO$Score)$p.value
  #Integrate conservation data
  cancer_conservation_data <- rbind(cancer_phastCons100way_vertebrate_box_data)
  cancer_mutation_method_wilcox_P1 <- data.frame(c(cancer_mutation_phastCons100way_vertebrate_L_LO_P,
                                                   cancer_mutation_phastCons100way_vertebrate_R_RO_P))
  colnames(cancer_mutation_method_wilcox_P1) <- c("P")
  cancer_mutation_method_wilcox_P$Method <- c('phastCons100way_vertebrate_L_LO',
                                              'phastCons100way_vertebrate_R_RO')
  cancer_mutation_method_wilcox_P1$p_lable <- ifelse(cancer_mutation_method_wilcox_P1$P>0.05,"NS",
                                                     ifelse(cancer_mutation_method_wilcox_P1$P>0.01,"*",
                                                            ifelse(cancer_mutation_method_wilcox_P1$P>0.001,"**",
                                                                   ifelse(cancer_mutation_method_wilcox_P1$P>0.0001,"***","****"))))
  if(plot1==TRUE){
    A_data <- data.frame(c(unique(cancer_damage_data$Method)))
    colnames(A_data) <- "Method"
    A_data$"Hugo_Symbol" <- 'A'
    A_data$"Score" <- 0
    A_data$"lable" <- 'A'
    A_data <- A_data[,colnames(cancer_damage_data)]
    cancer_damage_data <- rbind(cancer_damage_data,A_data)
    cancer_damage_data$lable<-factor(cancer_damage_data$lable,levels = c("ligand","ligand_other","A","receptor","receptor_other"))
    cancer_damage_data$Method <- factor(cancer_damage_data$Method,levels = c('Polyphen2 HDIV','Polyphen2 HVAR','DEOGEN2','LRT','M CAP','MetaLR'))
    p1 <- ggplot(cancer_damage_data, aes(x = Method, y = Score))+
      geom_boxplot(aes(fill=lable,color=lable),
                   alpha = 0.6# 透明度
                   ,lwd=1,width=0.8,
                   position = position_dodge(0.9),
                   outlier.shape = NA# 外点颜色
      )+
      scale_fill_manual(values = c(ligand="#1bb2b1",A=alpha("white",0),ligand_other="#bfbebe",receptor="#f4b248",receptor_other="#bfbebe"),
                        labels=c("ligand","ligand_other","","receptor","receptor_other"))+
      scale_color_manual(values = c(ligand="#1bb2b1",A=alpha("white",0),ligand_other="#bfbebe",receptor="#f4b248",receptor_other="#bfbebe"),
                         labels=c("ligand","ligand_other","","receptor","receptor_other"))+
      scale_y_continuous(breaks=seq(0,ylim1,ylim1/6))+
      coord_cartesian(ylim =  c(0,ylim1))+
      labs(x="Method",y="Prediction score")+  # 白色主题+
      geom_signif(annotations = cancer_mutation_method_wilcox_P$p_lable,
                  y_position = y_position1,
                  xmin = c(0.65,1.2,1.65,2.2,2.65,3.2,3.65,4.2,4.65,5.2,5.65,6.2),
                  xmax = c(0.80,1.35,1.8,2.35,2.8,3.35,3.8,4.35,4.8,5.35,5.8,6.35),
                  textsize=5,size = 0.8,tip_length = 0.02)+
      theme(axis.title = element_text(face = "bold",size = 20), # 坐标轴标题大小
            axis.text.y  = element_text(face = "bold",size = 12), # 坐标轴标签大小
            axis.text.x  = element_text(face = "bold",size = 12), # 坐标轴标签大小
            plot.title = element_text(size = 20,hjust = 0.5,face = "bold"), # 标题设置
            panel.background = element_blank(),#去除背景色和网格线
            panel.border = element_rect(fill = "NA",color="black", size =2),
            legend.position = "top",
            legend.direction = "horizontal"
      )
    print(p1)
  }
  if(plot2==TRUE){
    B_data <- data.frame(c(unique(cancer_conservation_data$Method)))
    colnames(B_data) <- "Method"
    B_data$"Hugo_Symbol" <- 'B'
    B_data$"Score" <- 0
    B_data$"lable" <- 'B'
    B_data <- B_data[,colnames(cancer_conservation_data)]
    cancer_conservation_data <- rbind(cancer_conservation_data,B_data)
    cancer_conservation_data$lable<-factor(cancer_conservation_data$lable,levels = c("ligand","ligand_other","B","receptor","receptor_other"))
    p2 <- ggplot(cancer_conservation_data, aes(x = Method, y = Score))+
      geom_boxplot(aes(fill=lable,color=lable),
                   alpha = 0.6
                   ,lwd=1,width=0.8,
                   position = position_dodge(0.9),
                   outlier.shape = NA
      )+
      scale_fill_manual(values = c(ligand="#1bb2b1",B=alpha("white",0),ligand_other="#bfbebe",receptor="#f4b248",receptor_other="#bfbebe"),
                        labels=c("ligand","ligand_other","","receptor","receptor_other"))+
      scale_color_manual(values = c(ligand="#1bb2b1",B=alpha("white",0),ligand_other="#bfbebe",receptor="#f4b248",receptor_other="#bfbebe"),
                         labels=c("ligand","ligand_other","","receptor","receptor_other"))+
      scale_y_continuous(breaks=seq(0,ylim2,ylim2/6))+
      coord_cartesian(ylim =  c(0,ylim2))+
      labs(x="phastCons100way vertebrated",y='Conservation score')+
      geom_signif(annotations = cancer_mutation_method_wilcox_P1$p_lable,
                  y_position = y_position2,
                  xmin = c(0.65,1.2),
                  xmax = c(0.80,1.35),
                  textsize=5,size = 0.8,tip_length = 0.02)+
      theme(axis.title = element_text(face = "bold",size = 20),
            axis.text.x = element_blank(),
            axis.text.y  = element_text(face = "bold",size = 12),
            plot.title = element_text(size = 20,hjust = 0.5,face = "bold"),
            panel.background = element_blank(),
            panel.border = element_rect(fill = "NA",color="black", size =2),
            legend.position = "top",
            legend.direction = "horizontal"
      )
    print(p2)
  }
  if(plot1==F&plot2==F)
    return(list(cancer_damage_data,cancer_mutation_method_wilcox_P,cancer_conservation_data,cancer_mutation_method_wilcox_P1))
}
