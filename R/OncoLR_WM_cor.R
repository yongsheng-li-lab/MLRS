#' @title Characterizing the association of oncoLRI in mutation and wild-type.
#'
#' @param OncoLRI Preferred oncoLRI.
#' @param exp Cancer Gene Expression Profiling File, FPKM Values.
#' @param mutation_data A dataframe, defaults to TCGA-MC3 project data.
#' @param ligand.source.cell Ligand source cells. Default 'Tumor cell'
#' @param receptor.source.cell Receptor source cells. Default 'Other cell'
#'
#' @return
#' @export
#'
#' @examples
OncoLR_WM_cor <- function(OncoLR,exp,mutation_data,ligand.source.cell='Tumor cell',receptor.source.cell='Other cell'){
  ######Reading mutation data
  print('Step1-------Calculating mutation')
  gene <- as.character(unique(mutation_data$Hugo_Symbol))
  sample <- as.character(unique(mutation_data$Tumor_Sample_Barcode))
  matutation_0_1 <- as.data.frame(matrix(0,length(gene),length(sample),
                                              dimnames = list(gene,sample)))
  for (i in 1:nrow(mutation_data)){
    matutation_0_1[as.character(mutation_data[i,"Hugo_Symbol"]),as.character(mutation_data[i,"Tumor_Sample_Barcode"])] <- 1
  }
  ######Reading expression profile
  ###Calculating mutation
  print('Step2-------Reading expression profile')
  colnames(exp) <- substr(colnames(exp),1,16)
  sample_id <- colnames(exp)
  tumor_sample_ID=grep(".*-0\\d[ABC]",sample_id,value = T)
  colnames(matutation_0_1) <- substr(colnames(matutation_0_1),1,16)
  exp_mutation_sample <- intersect(colnames(exp),colnames(matutation_0_1))#1022个样本
  exp_mutation_tumor_sample <- intersect(exp_mutation_sample,tumor_sample_ID)
  exp1 <- exp[,c(exp_mutation_tumor_sample)]
  matutation_0_1_1 <- matutation_0_1[,c(exp_mutation_tumor_sample)]
  matutation_0_1_2 <- data.frame(t(matutation_0_1_1))
  colnames(matutation_0_1_2) <- gsub(r'(\.)',"-",colnames(matutation_0_1_2))
  ######Run OncoLR_WM_cor
  print('Step3-------Running OncoLR_WM_cor')
  exp2 <- log2(exp1+1)
  if(ligand.source.cell=='Tumor cell'){
    LWR_cor_res <- data.frame()
    for(i in 1:nrow(OncoLR)){
      wild_sample_L <- rownames(matutation_0_1_2[which(matutation_0_1_2[,OncoLR[i,'ligand']]==0),])
      if(length(wild_sample_L)<=5){
        cor_p=NA
        p_p=NA
        cor_s=NA
        p_s=NA
        x=c(OncoLR[i,'ligand'],
            OncoLR[i,'receptor'],
            cor_p,p_p,cor_s,p_s)
        LWR_cor_res=rbind(LWR_cor_res,x)
        cor_p=c()
        p_p=c()
        cor_s=c()
        p_s=c()
      }else{
        a=exp2[OncoLR[i,'ligand'],wild_sample_L]
        b=exp2[OncoLR[i,'receptor'],wild_sample_L]
        cor_p=cor.test(as.numeric(a),as.numeric(b),method = "pearson")[["estimate"]][["cor"]]
        p_p=cor.test(as.numeric(a),as.numeric(b),method = "pearson")[["p.value"]]
        cor_s=cor.test(as.numeric(a),as.numeric(b),method = "spearman")[["estimate"]][["rho"]]
        p_s=cor.test(as.numeric(a),as.numeric(b),method = "spearman")[["p.value"]]
        x=c(OncoLR[i,'ligand'],
            OncoLR[i,'receptor'],
            cor_p,p_p,cor_s,p_s)
        LWR_cor_res=rbind(LWR_cor_res,x)
        a=c()
        b=c()
        cor_p=c()
        p_p=c()
        cor_s=c()
        p_s=c()
      }
    }
    colnames(LWR_cor_res) <- c('ligand','receptor','pearson_cor','pearson_p','spearman_cor','spearman_p')
    LWR_cor_res <- na.omit(LWR_cor_res)
    LWR_cor_res$pearson_cor <- as.numeric(LWR_cor_res$pearson_cor)
    LWR_cor_res$pearson_p <- as.numeric(LWR_cor_res$pearson_p)
    LWR_cor_res$pearson_FDR <- p.adjust(LWR_cor_res$pearson_p,method = 'BH',n=length(LWR_cor_res$pearson_p))
    LWR_cor_res$spearman_cor <- as.numeric(LWR_cor_res$spearman_cor)
    LWR_cor_res$spearman_p <- as.numeric(LWR_cor_res$spearman_p)
    LWR_cor_res$spearman_FDR <- p.adjust(LWR_cor_res$spearman_p,method = 'BH',n=length(LWR_cor_res$spearman_p))
    LMR_cor_res <- data.frame()
    for(i in 1:nrow(OncoLR)){
      mutation_sample_L <- rownames(matutation_0_1_2[which(matutation_0_1_2[,OncoLR[i,'ligand']]==1),])
      n=length(mutation_sample_L)
      if(n<=5){
        cor_p=NA
        p_p=NA
        cor_s=NA
        p_s=NA
        x=c(OncoLR[i,'ligand'],
            OncoLR[i,'receptor'],
            cor_p,p_p,cor_s,p_s)
        LMR_cor_res=rbind(LMR_cor_res,x)
        cor_p=c()
        p_p=c()
        cor_s=c()
        p_s=c()
      }else{
        a=exp2[OncoLR[i,'ligand'],mutation_sample_L]
        b=exp2[OncoLR[i,'receptor'],mutation_sample_L]
        cor_p=cor.test(as.numeric(a),as.numeric(b),method = "pearson")[["estimate"]][["cor"]]
        p_p=cor.test(as.numeric(a),as.numeric(b),method = "pearson")[["p.value"]]
        cor_s=cor.test(as.numeric(a),as.numeric(b),method = "spearman")[["estimate"]][["rho"]]
        p_s=cor.test(as.numeric(a),as.numeric(b),method = "spearman")[["p.value"]]
        x=c(OncoLR[i,'ligand'],
            OncoLR[i,'receptor'],
            cor_p,p_p,cor_s,p_s)
        LMR_cor_res=rbind(LMR_cor_res,x)
        a=c()
        b=c()
        cor.r=c()
        p=c()
      }
    }
    colnames(LMR_cor_res) <- c('ligand','receptor','pearson_cor','pearson_p','spearman_cor','spearman_p')
    LMR_cor_res <- na.omit(LMR_cor_res)
    LMR_cor_res$pearson_cor <- as.numeric(LMR_cor_res$pearson_cor)
    LMR_cor_res$pearson_p <- as.numeric(LMR_cor_res$pearson_p)
    LMR_cor_res$pearson_FDR <- p.adjust(LMR_cor_res$pearson_p,method = 'BH',n=length(LMR_cor_res$pearson_p))
    LMR_cor_res$spearman_cor <- as.numeric(LMR_cor_res$spearman_cor)
    LMR_cor_res$spearman_p <- as.numeric(LMR_cor_res$spearman_p)
    LMR_cor_res$spearman_FDR <- p.adjust(LMR_cor_res$spearman_p,method = 'BH',n=length(LMR_cor_res$spearman_p))
    #####Return
    print('Step4-------OncoLR_WM_cor is over')
    LWR_cor_res1 <- LWR_cor_res
    colnames(LWR_cor_res1) <- c('ligand','receptor',"wild_pearson_cor","wild_pearson_p","wild_spearman_cor","wild_spearman_p","wild_pearson_FDR","wild_spearman_FDR")
    LWR_cor_res1$interaction <- paste0(LWR_cor_res1$ligand,"_",LWR_cor_res1$receptor)
    LMR_cor_res1 <- LMR_cor_res
    colnames(LMR_cor_res1) <- c('ligand','receptor',"mutation_pearson_cor","mutation_pearson_p","mutation_spearman_cor","mutation_spearman_p","mutation_pearson_FDR","mutation_spearman_FDR")
    LMR_cor_res1$interaction <- paste0(LMR_cor_res1$ligand,"_",LMR_cor_res1$receptor)
    TLR_cor_data1 <- merge(LWR_cor_res1,LMR_cor_res1,by='interaction')
    TLR_cor_data1 <- TLR_cor_data1[,c('ligand.x','receptor.x','interaction',"wild_pearson_cor","wild_pearson_p","wild_spearman_cor","wild_spearman_p","wild_pearson_FDR","wild_spearman_FDR","mutation_pearson_cor","mutation_pearson_p","mutation_spearman_cor","mutation_spearman_p","mutation_pearson_FDR","mutation_spearman_FDR")]
    colnames(TLR_cor_data1) <- c('ligand','receptor','interaction',"wild_pearson_cor","wild_pearson_p","wild_spearman_cor","wild_spearman_p","wild_pearson_FDR","wild_spearman_FDR","mutation_pearson_cor","mutation_pearson_p","mutation_spearman_cor","mutation_spearman_p","mutation_pearson_FDR","mutation_spearman_FDR")
    TLR_cor_data1 <- na.omit(TLR_cor_data1)
    TLR_cor_data1$P_lable <- ifelse(TLR_cor_data1$wild_pearson_p<0.05&TLR_cor_data1$mutation_pearson_p<0.05,'Sig-Sig',
                                         ifelse(TLR_cor_data1$wild_pearson_p<0.05&TLR_cor_data1$mutation_pearson_p>0.05,'Sig-No',
                                                ifelse(TLR_cor_data1$wild_pearson_p>0.05&TLR_cor_data1$mutation_pearson_p<0.05,'No-Sig','No-No')))
    return(TLR_cor_data1)
  }
  if(receptor.source.cell=='Tumor cell'){
    LRW_cor_res <- data.frame()
    for(i in 1:nrow(OncoLR)){
      wild_sample_R <- rownames(matutation_0_1_2[which(matutation_0_1_2[,OncoLR[i,'receptor']]==0),])
      if(length(wild_sample_R)<=5){
        cor_p=NA
        p_p=NA
        cor_s=NA
        p_s=NA
        x=c(OncoLR[i,'ligand'],
            OncoLR[i,'receptor'],
            cor_p,p_p,cor_s,p_s)
        LRW_cor_res=rbind(LRW_cor_res,x)
        cor_p=c()
        p_p=c()
        cor_s=c()
        p_s=c()
      }else{
        a=exp2[OncoLR[i,'ligand'],wild_sample_R]
        b=exp2[OncoLR[i,'receptor'],wild_sample_R]
        cor_p=cor.test(as.numeric(a),as.numeric(b),method = "pearson")[["estimate"]][["cor"]]
        p_p=cor.test(as.numeric(a),as.numeric(b),method = "pearson")[["p.value"]]
        cor_s=cor.test(as.numeric(a),as.numeric(b),method = "spearman")[["estimate"]][["rho"]]
        p_s=cor.test(as.numeric(a),as.numeric(b),method = "spearman")[["p.value"]]
        x=c(OncoLR[i,'ligand'],
            OncoLR[i,'receptor'],
            cor_p,p_p,cor_s,p_s)
        LRW_cor_res=rbind(LRW_cor_res,x)
        a=c()
        b=c()
        cor_p=c()
        p_p=c()
        cor_s=c()
        p_s=c()
      }
    }
    colnames(LRW_cor_res) <- c('ligand','receptor','pearson_cor','pearson_p','spearman_cor','spearman_p')
    LRW_cor_res <- na.omit(LRW_cor_res)
    LRW_cor_res$pearson_cor <- as.numeric(LRW_cor_res$pearson_cor)
    LRW_cor_res$pearson_p <- as.numeric(LRW_cor_res$pearson_p)
    LRW_cor_res$pearson_FDR <- p.adjust(LRW_cor_res$pearson_p,method = 'BH',n=length(LRW_cor_res$pearson_p))
    LRW_cor_res$spearman_cor <- as.numeric(LRW_cor_res$spearman_cor)
    LRW_cor_res$spearman_p <- as.numeric(LRW_cor_res$spearman_p)
    LRW_cor_res$spearman_FDR <- p.adjust(LRW_cor_res$spearman_p,method = 'BH',n=length(LRW_cor_res$spearman_p))
    LRM_cor_res <- data.frame()
    for(i in 1:nrow(OncoLR)){
      mutation_sample_R <- rownames(matutation_0_1_2[which(matutation_0_1_2[,OncoLR[i,'receptor']]==1),])
      n=length(mutation_sample_R)
      if(n<=5){
        cor_p=NA
        p_p=NA
        cor_s=NA
        p_s=NA
        x=c(OncoLR[i,'ligand'],
            OncoLR[i,'receptor'],
            cor_p,p_p,cor_s,p_s)
        LRM_cor_res=rbind(LRM_cor_res,x)
        cor_p=c()
        p_p=c()
        cor_s=c()
        p_s=c()
      }else{
        a=exp2[OncoLR[i,'ligand'],mutation_sample_R]
        b=exp2[OncoLR[i,'receptor'],mutation_sample_R]
        cor_p=cor.test(as.numeric(a),as.numeric(b),method = "pearson")[["estimate"]][["cor"]]
        p_p=cor.test(as.numeric(a),as.numeric(b),method = "pearson")[["p.value"]]
        cor_s=cor.test(as.numeric(a),as.numeric(b),method = "spearman")[["estimate"]][["rho"]]
        p_s=cor.test(as.numeric(a),as.numeric(b),method = "spearman")[["p.value"]]
        x=c(OncoLR[i,'ligand'],
            OncoLR[i,'receptor'],
            cor_p,p_p,cor_s,p_s)
        LRM_cor_res=rbind(LRM_cor_res,x)
        a=c()
        b=c()
        cor.r=c()
        p=c()
      }
    }
    colnames(LRM_cor_res) <- c('ligand','receptor','pearson_cor','pearson_p','spearman_cor','spearman_p')
    LRM_cor_res <- na.omit(LRM_cor_res)
    LRM_cor_res$pearson_cor <- as.numeric(LRM_cor_res$pearson_cor)
    LRM_cor_res$pearson_p <- as.numeric(LRM_cor_res$pearson_p)
    LRM_cor_res$pearson_FDR <- p.adjust(LRM_cor_res$pearson_p,method = 'BH',n=length(LRM_cor_res$pearson_p))
    LRM_cor_res$spearman_cor <- as.numeric(LRM_cor_res$spearman_cor)
    LRM_cor_res$spearman_p <- as.numeric(LRM_cor_res$spearman_p)
    LRM_cor_res$spearman_FDR <- p.adjust(LRM_cor_res$spearman_p,method = 'BH',n=length(LRM_cor_res$spearman_p))
    LRT_cor_data <- rbind(LRW_cor_res,LRM_cor_res)
    LRT_cor_data$Mutation <- c(rep('Wild',nrow(LRW_cor_res)),
                               rep('Mutation',nrow(LRM_cor_res)))
    #####Return
    print('Step4-------OncoLR_WM_cor is over')
    LRW_cor_res1  <- LRW_cor_res
    colnames(LRW_cor_res1 ) <- c('ligand','receptor',"wild_pearson_cor","wild_pearson_p","wild_spearman_cor","wild_spearman_p","wild_pearson_FDR","wild_spearman_FDR")
    LRW_cor_res1$interaction <- paste0(LRW_cor_res1 $ligand,"_",LRW_cor_res1 $receptor)
    LRM_cor_res1 <- LRM_cor_res
    colnames(LRM_cor_res1) <- c('ligand','receptor',"mutation_pearson_cor","mutation_pearson_p","mutation_spearman_cor","mutation_spearman_p","mutation_pearson_FDR","mutation_spearman_FDR")
    LRM_cor_res1$interaction <- paste0(LRM_cor_res1$ligand,"_",LRM_cor_res1$receptor)
    LRT_cor_data1 <- merge(LRW_cor_res1 ,LRM_cor_res1,by='interaction')
    LRT_cor_data1 <- LRT_cor_data1[,c('ligand.x','receptor.x','interaction',"wild_pearson_cor","wild_pearson_p","wild_spearman_cor","wild_spearman_p","wild_pearson_FDR","wild_spearman_FDR","mutation_pearson_cor","mutation_pearson_p","mutation_spearman_cor","mutation_spearman_p","mutation_pearson_FDR","mutation_spearman_FDR")]
    colnames(LRT_cor_data1) <- c('ligand','receptor','interaction',"wild_pearson_cor","wild_pearson_p","wild_spearman_cor","wild_spearman_p","wild_pearson_FDR","wild_spearman_FDR","mutation_pearson_cor","mutation_pearson_p","mutation_spearman_cor","mutation_spearman_p","mutation_pearson_FDR","mutation_spearman_FDR")
    LRT_cor_data1 <- na.omit(LRT_cor_data1)
    LRT_cor_data1$P_lable <- ifelse(LRT_cor_data1$wild_pearson_p<0.05&LRT_cor_data1$mutation_pearson_p<0.05,'Sig-Sig',
                                         ifelse(LRT_cor_data1$wild_pearson_p<0.05&LRT_cor_data1$mutation_pearson_p>0.05,'Sig-No',
                                                ifelse(LRT_cor_data1$wild_pearson_p>0.05&LRT_cor_data1$mutation_pearson_p<0.05,'No-Sig','No-No')))
    return(LRT_cor_data1)
    }
}
