#' @title Preferred ligand-receptor pairs.
#'
#' @param CCC_net Cell-cell communication network.
#' @param exp1 Cancer Gene Expression Profiling File 1, FPKM Values.
#' @param exp2 Cancer gene expression profile file 2, FPKM_UQ values.
#' @param mutation_data a dataframe, defaults to TCGA-MC3 project data.
#' @param q The preference threshold, which determines whether the ligand-receptor pair is significant, has a default value of 0.01.
#'
#' @importFrom aod
#' @importFrom estimate
#' @importFrom fdrtool
#'
#' @return
#' @export
#'
#' @examples
getOncoLRI <- function(CCC_net,exp1,exp2,mutation_data,q=0.01){
  ######Reading mutation data
  ###Calculating mutation
  print('Step1-------Calculating mutation')
  gene <- as.character(unique(mutation_data$Hugo_Symbol))
  sample <- as.character(unique(mutation_data$Tumor_Sample_Barcode))
  matutation_0_1 <- as.data.frame(matrix(0,length(gene),length(sample),
                                              dimnames = list(gene,sample)))
  for (i in 1:nrow(mutation_data)){
    matutation_0_1[as.character(mutation_data[i,"Hugo_Symbol"]),as.character(mutation_data[i,"Tumor_Sample_Barcode"])] <- 1
  }
  ######Reading expression profile
  print('Step2-------Reading expression profile')
  sample_id <- colnames(exp1)
  exp_mutation_tumor_sample <- intersect(colnames(exp1),colnames(matutation_0_1))
  exp11 <- exp1[,c(exp_mutation_tumor_sample)]
  matutation_0_1_1 <- matutation_0_1[,c(exp_mutation_tumor_sample)]
  #####Reading CCC net
  print('Step3-------Reading CCC net')
  ##Separate out the interactions between cancer cells and other cells, and between other cells and cancer cells
  ligand <- data.frame(intersect(unique(CCC_net$ligand),rownames(exp11)))
  colnames(ligand) <- c('ligand')
  receptor <- data.frame(intersect(unique(CCC_net$receptor),rownames(exp11)))
  colnames(receptor) <- c('receptor')
  CCC_net_exp1 <- merge(CCC_net,ligand,by='ligand')
  CCC_net_exp2 <- merge(CCC_net_exp1,receptor,by='receptor')
  CCC_net_exp2 <- CCC_net_exp2[!duplicated(CCC_net_exp2),]
  CCC_net_exp <- CCC_net_exp2[,c("source","ligand","receptor","target","interaction_name")]
  TCC_net_exp <- CCC_net_exp[which(CCC_net_exp$source=='Tumor'),]
  CTC_net_exp <- CCC_net_exp[which(CCC_net_exp$target=='Tumor'),]
  ######Repurposing cancer mutation data
  print('Step4-------Repurposing cancer mutation data')
  matutation_0_1_2 <- data.frame(t(matutation_0_1_1))
  colnames(matutation_0_1_2) <- gsub(r'(\.)',"-",colnames(matutation_0_1_2))
  ######Extraction of L-R interactions with mutational information in tumor cell-derived L-R interactions
  LT_mutation <- data.frame(intersect(unique(TCC_net_exp$ligand),colnames(matutation_0_1_2)))
  colnames(LT_mutation) <- c('ligand')
  LTCC_LR_mutation <- merge(LT_mutation,TCC_net_exp,by='ligand')
  LTCC_LR_mutation <- LTCC_LR_mutation[!duplicated(LTCC_LR_mutation),]
  LTCC_LR_mutation <- LTCC_LR_mutation[,c("source","ligand","receptor","target","interaction_name")]
  #####Extraction of interactions with mutational information in L-R tumor cell-derived interactions
  RT_mutation <- data.frame(intersect(unique(CTC_net_exp$receptor),colnames(matutation_0_1_2)))
  colnames(RT_mutation) <- c('receptor')
  CRTC_LR_mutation <- merge(RT_mutation,CTC_net_exp,by='receptor')
  CRTC_LR_mutation <- CRTC_LR_mutation[!duplicated(CRTC_LR_mutation),]
  CRTC_LR_mutation <- CRTC_LR_mutation[,c("source","ligand","receptor","target","interaction_name")]
  ##########Preferences were made using a regression model, this time preferring to carry cell labels
  print('Step5-------Calculate substrate score')
  #######Calculation of the substrate score was first carried out
  ###Calculate substrate score-remember to convert to FPKM-UQ values
  ##Organizing Expression
  exp21 <- exp2[,exp_mutation_tumor_sample]
  exp21$id <- rownames(exp21)
  exp21 <- exp21[,c(ncol(exp21),1:(ncol(exp21)-1))]
  write.table(exp21,'exp21.txt',row.names = F,sep = '\t',quote = F)
  ##Calculate substrate score
  filterCommonGenes(input.f="exp21.txt", output.f="exp21.gct", id="GeneSymbol")
  estimateScore(input.ds="exp21.gct", output.ds="estimate_result.gct", platform="affymetrix")
  estimate_result <- read.table("estimate_result.gct", skip = 2, header = TRUE)
  rownames(estimate_result) <- estimate_result$NAME
  estimate_result<- estimate_result[,-c(1,2)]
  TumorPurity_score <- data.frame(estimate_result['TumorPurity',])
  colnames(TumorPurity_score) <- gsub(r'(\.)',"-",colnames(TumorPurity_score))
  TumorPurity_score <- TumorPurity_score[,colnames(exp11)]
  Tumor_Purity <- data.frame(t(TumorPurity_score))
  #######Introducing a linear model with an interaction term that determines whether the correlation strengthens or weakens based on the positivity or negativity of β1, and a q-value that determines whether the event is significant or not
  ####Write a loop to count the tumor cell-LR pairs affected by mutations between them
  ##At this point we consider 2 scenarios
  #L/R would appear to have an expression profile of 0, so on this line we go all the way to NA
  #β3 will appear NA, this is because the variables appear strongly correlated, i.e. cor=1,so in this row we all become NA
  ##Run linear model
  print('Step6-------Run linear model')
  exp12 <- log2(exp11+1)
  ##Computing L(tumor cell)-R derived interactions to extract interactions with mutational information
  TLR_model_res <- data.frame()
  for(i in 1:nrow(LTCC_LR_mutation)){
    exp_data <- data.frame(t(exp12[c(LTCC_LR_mutation[i,"ligand"],LTCC_LR_mutation[i,"receptor"]),]))
    if(sum(exp_data[,1])==0||sum(exp_data[,2])==0){
      b1 <- NA
      b2 <- NA
      b3 <- NA
      p2 <- NA
      p3 <- NA
      p_value <- NA
      x <- data.frame(LTCC_LR_mutation[i,],b1,b2,b3,p2,p3,p_value)
      TLR_model_res <- rbind(TLR_model_res,x)
      exp_data <- data.frame()
      Mutation_Status <- data.frame()
      data <- data.frame()
      x <- data.frame()
    }else{
      Mutation_Status <- data.frame(t(matutation_0_1_1[LTCC_LR_mutation[i,"ligand"],]))
      ligandMutation <- exp_data[,1]*Mutation_Status
      data <- cbind(exp_data,Mutation_Status,ligandMutation,Tumor_Purity)
      colnames(data) <- c('ligand','receptor','Mutation_Status','ligandMutation','Tumor_Purity')
      model <- lm(receptor ~ ligand+Mutation_Status+ligandMutation+Tumor_Purity, data = data)
      b1 <- model[["coefficients"]][["ligand"]]
      b2 <- model[["coefficients"]][["Mutation_Status"]]
      b3 <- model[["coefficients"]][["ligandMutation"]]
      if(is.na(b2)==TRUE||is.na(b3)==TRUE){
        p2 <- NA
        p3 <- NA
        p_value <- NA
        x <- data.frame(LTCC_LR_mutation[i,],b1,b2,b3,p2,p3,p_value)
        TLR_model_res <- rbind(TLR_model_res,x)
        exp_data <- data.frame()
        Mutation_Status <- data.frame()
        data <- data.frame()
        x <- data.frame()
      }else{
        p2 <- wald.test(b = model[["coefficients"]], Sigma = vcov(model), Terms = 2)[["result"]][["chi2"]][["P"]]
        p3 <- wald.test(b = model[["coefficients"]], Sigma = vcov(model), Terms = 3)[["result"]][["chi2"]][["P"]]
        p_values <- c(p2, p3)
        log_p_values <- log(p_values)
        sum_log_p <- -2*sum(log_p_values)
        p_value <- 1 - pchisq(sum_log_p, df = 2 * length(p_values))
        x <- data.frame(LTCC_LR_mutation[i,],b1,b2,b3,p2,p3,p_value)
        TLR_model_res <- rbind(TLR_model_res,x)
        exp_data <- data.frame()
        Mutation_Status <- data.frame()
        data <- data.frame()
        x <- data.frame()
      }
    }
  }
  TLR_model_res1 <- na.omit(TLR_model_res)
  TLR_model_res1$q <- fdrtool(TLR_model_res1$p_value,statistic = "pvalue")$qval
  TLR_model_res1$Condition <- ifelse(TLR_model_res1$q<q,'Sig','No-sig')
  TLR_model_res2 <- TLR_model_res1[which(TLR_model_res1$q<q),]
  ##Computing L-R tumor cell-derived interactions to extract interactions with mutational information
  LRT_model_res <- data.frame()
  for(i in 1:nrow(CRTC_LR_mutation)){
    exp_data <- data.frame(t(exp12[c(CRTC_LR_mutation[i,"ligand"],CRTC_LR_mutation[i,"receptor"]),]))
    if(sum(exp_data[,1])==0||sum(exp_data[,2])==0){
      b1 <- NA
      b2 <- NA
      b3 <- NA
      p2 <- NA
      p3 <- NA
      p_value <- NA
      x <- data.frame(CRTC_LR_mutation[i,],b1,b2,b3,p2,p3,p_value)
      LRT_model_res <- rbind(LRT_model_res,x)
      exp_data <- data.frame()
      Mutation_Status <- data.frame()
      data <- data.frame()
      x <- data.frame()
    }else{
      Mutation_Status <- data.frame(t(matutation_0_1_1[CRTC_LR_mutation[i,"receptor"],]))
      receptorMutation <- exp_data[,2]*Mutation_Status
      data <- cbind(exp_data,Mutation_Status,receptorMutation,Tumor_Purity)
      colnames(data) <- c('ligand','receptor','Mutation_Status','receptorMutation','Tumor_Purity')
      model <- lm(ligand ~ receptor +Mutation_Status+receptorMutation+Tumor_Purity, data = data)
      b1 <- model[["coefficients"]][["receptor"]]
      b2 <- model[["coefficients"]][["Mutation_Status"]]
      b3 <- model[["coefficients"]][["receptorMutation"]]
      if(is.na(b2)==TRUE||is.na(b3)==TRUE){
        p2 <- NA
        p3 <- NA
        p_value <- NA
        x <- data.frame(CRTC_LR_mutation[i,],b1,b2,b3,p2,p3,p_value)
        LRT_model_res <- rbind(LRT_model_res,x)
        exp_data <- data.frame()
        Mutation_Status <- data.frame()
        data <- data.frame()
        x <- data.frame()
      }else{
        p2 <- wald.test(b = model[["coefficients"]], Sigma = vcov(model), Terms = 2)[["result"]][["chi2"]][["P"]]
        p3 <- wald.test(b = model[["coefficients"]], Sigma = vcov(model), Terms = 3)[["result"]][["chi2"]][["P"]]
        p_values <- c(p2, p3)
        log_p_values <- log(p_values)
        sum_log_p <- -2*sum(log_p_values)
        p_value <- 1 - pchisq(sum_log_p, df = 2 * length(p_values))
        x <- data.frame(CRTC_LR_mutation[i,],b1,b2,b3,p2,p3,p_value)
        LRT_model_res <- rbind(LRT_model_res,x)
        exp_data <- data.frame()
        Mutation_Status <- data.frame()
        data <- data.frame()
        x <- data.frame()
      }
    }
  }
  LRT_model_res1 <- na.omit(LRT_model_res)
  LRT_model_res1$q <- fdrtool(LRT_model_res1$p_value,statistic = "pvalue")$qval
  LRT_model_res1$Condition <- ifelse(LRT_model_res1$q<q,'Sig','No-sig')
  LRT_model_res2 <- LRT_model_res1[which(LRT_model_res1$q<q),]
  ############Return
  print('Step7-------Linear model is over')
  TLR_model_result <- TLR_model_res2
  LRT_model_result <- LRT_model_res2
  OncoLRI <- list(TLR_model_result,LRT_model_result)
  return(OncoLRI)
}
