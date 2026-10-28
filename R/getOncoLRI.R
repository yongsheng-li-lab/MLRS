#' @title Preferred ligand-receptor pairs.
#'
#' @param CCC_net Cell-cell communication network.
#' @param exp1 Cancer Gene Expression Profiling File 1, FPKM Values.
#' @param exp2 Cancer gene expression profile file 2, FPKM_UQ Values.
#' @param mutation_data Mutation data containing Hugo_Symbol and Tumor_Sample_Barcode.
#' @param q Significance threshold for the q value, default 0.01.
#'
#' @return A list with TLR and LRT significant cell-resolved interactions.
#' @export
#'
#' @examples
#' OncoLRI <- getOncoLRI(CCC_net, exp1, exp2, mutation_data, q = 0.01)
getOncoLRI <- function(CCC_net,exp1,exp2,mutation_data,q=0.01){
  ######Reading mutation data
  ###Calculating mutation
  print('Step1-------Calculating mutation')
  gene <- as.character(unique(mutation_data$Hugo_Symbol))
  sample <- as.character(unique(mutation_data$Tumor_Sample_Barcode))
  matutation_0_1 <- as.data.frame(
    matrix(
      0,
      length(gene),
      length(sample),
      dimnames = list(gene,sample)
    )
  )
  for(i in 1:nrow(mutation_data)){
    matutation_0_1[
      as.character(mutation_data[i,"Hugo_Symbol"]),
      as.character(mutation_data[i,"Tumor_Sample_Barcode"])
    ] <- 1
  }
  ######Reading expression profile
  print('Step2-------Reading expression profile')
  ###Use only exact sample IDs shared by FPKM, FPKM-UQ and mutation data.
  ###No automatic TCGA/CPTAC ID conversion is performed here.
  exp_mutation_tumor_sample <- Reduce(
    intersect,
    list(
      colnames(exp1),
      colnames(exp2),
      colnames(matutation_0_1)
    )
  )
  if(length(exp_mutation_tumor_sample) < 10){
    stop(
      paste0(
        "Fewer than 10 exact matched samples among exp1, exp2 and mutation_data: ",
        length(exp_mutation_tumor_sample)
      )
    )
  }
  exp11 <- exp1[,exp_mutation_tumor_sample,drop=FALSE]
  exp21 <- exp2[,exp_mutation_tumor_sample,drop=FALSE]
  matutation_0_1_1 <- matutation_0_1[,exp_mutation_tumor_sample,drop=FALSE]
  cat("Matched samples:",length(exp_mutation_tumor_sample),"\n")
  #####Reading CCC net
  print('Step3-------Reading CCC net')
  if(!all(c("source","ligand","receptor","target") %in% colnames(CCC_net))){
    stop("CCC_net must contain source, ligand, receptor and target")
  }
  if(!"interaction_name" %in% colnames(CCC_net)){
    CCC_net$interaction_name <- paste0(CCC_net$ligand,"_",CCC_net$receptor)
  }
  ###Separate interactions involving tumor cells
  ligand <- data.frame(
    ligand = intersect(unique(CCC_net$ligand),rownames(exp11)),
    stringsAsFactors = FALSE
  )
  receptor <- data.frame(
    receptor = intersect(unique(CCC_net$receptor),rownames(exp11)),
    stringsAsFactors = FALSE
  )
  CCC_net_exp1 <- merge(CCC_net,ligand,by="ligand")
  CCC_net_exp2 <- merge(CCC_net_exp1,receptor,by="receptor")
  CCC_net_exp2 <- CCC_net_exp2[!duplicated(CCC_net_exp2),,drop=FALSE]
  CCC_net_exp <- CCC_net_exp2[
    ,
    c("source","ligand","receptor","target","interaction_name"),
    drop=FALSE
  ]
  TCC_net_exp <- CCC_net_exp[CCC_net_exp$source=="Tumor",,drop=FALSE]
  CTC_net_exp <- CCC_net_exp[CCC_net_exp$target=="Tumor",,drop=FALSE]
  ######Repurposing cancer mutation data
  print('Step4-------Repurposing cancer mutation data')
  ###The mutation matrix already uses gene symbols as row names.
  ###Use them directly; no extra gene/sample name conversion is performed.
  mutation_genes <- rownames(matutation_0_1_1)
  #####Tumor-source direction: mutation is on ligand
  TLR_unique_LR <- unique(
    TCC_net_exp[,c("ligand","receptor"),drop=FALSE]
  )
  TLR_unique_LR <- TLR_unique_LR[
    TLR_unique_LR$ligand %in% mutation_genes,
    ,
    drop=FALSE
  ]
  #####Tumor-target direction: mutation is on receptor
  LRT_unique_LR <- unique(
    CTC_net_exp[,c("ligand","receptor"),drop=FALSE]
  )

  LRT_unique_LR <- LRT_unique_LR[
    LRT_unique_LR$receptor %in% mutation_genes,
    ,
    drop=FALSE
  ]

  cat("Tumor-source unique LR for regression:",nrow(TLR_unique_LR),"\n")
  cat("Tumor-target unique LR for regression:",nrow(LRT_unique_LR),"\n")

  ##########Calculate tumor purity
  print('Step5-------Calculate substrate score')

  exp21_estimate <- exp21
  exp21_estimate$id <- rownames(exp21_estimate)
  exp21_estimate <- exp21_estimate[
    ,
    c(ncol(exp21_estimate),1:(ncol(exp21_estimate)-1)),
    drop=FALSE
  ]

  ###Use a temporary directory so different runs do not overwrite each other
  estimate_tmp <- tempfile("getOncoLRI_estimate_")
  dir.create(estimate_tmp)
  on.exit(unlink(estimate_tmp,recursive=TRUE,force=TRUE),add=TRUE)

  exp21_txt <- file.path(estimate_tmp,"exp21.txt")
  exp21_gct <- file.path(estimate_tmp,"exp21.gct")
  estimate_gct <- file.path(estimate_tmp,"estimate_result.gct")

  write.table(
    exp21_estimate,
    exp21_txt,
    row.names=FALSE,
    sep="\t",
    quote=FALSE
  )

  estimate::filterCommonGenes(
    input.f=exp21_txt,
    output.f=exp21_gct,
    id="GeneSymbol"
  )

  estimate::estimateScore(
    input.ds=exp21_gct,
    output.ds=estimate_gct,
    platform="affymetrix"
  )

  estimate_result <- read.table(
    estimate_gct,
    skip=2,
    header=TRUE,
    check.names=FALSE
  )

  rownames(estimate_result) <- estimate_result$NAME
  estimate_result <- estimate_result[,-c(1,2),drop=FALSE]

  ###ESTIMATE/read.table can convert '-' to '.' in some R settings.
  ###Only restore that formatting; do not perform additional ID matching.
  colnames(estimate_result) <- gsub(r'(\.)',"-",colnames(estimate_result))

  if(!"TumorPurity" %in% rownames(estimate_result)){
    stop("TumorPurity was not found in the ESTIMATE result")
  }

  if(!all(colnames(exp11) %in% colnames(estimate_result))){
    stop("ESTIMATE sample IDs do not match the exact expression sample IDs")
  }

  TumorPurity_score <- data.frame(
    estimate_result["TumorPurity",colnames(exp11),drop=FALSE],
    check.names=FALSE
  )

  Tumor_Purity <- data.frame(t(TumorPurity_score))
  colnames(Tumor_Purity) <- "Tumor_Purity"

  ######Run linear model
  print('Step6-------Run linear model for unique LR pairs')

  exp12 <- log2(exp11+1)

  ########################################################################
  ##### TLR: Tumor is source; ligand mutation
  ########################################################################
  TLR_model_res <- data.frame()

  if(nrow(TLR_unique_LR)>0){
    for(i in seq_len(nrow(TLR_unique_LR))){

      LR_i <- TLR_unique_LR[i,,drop=FALSE]

      exp_data <- data.frame(
        t(exp12[c(LR_i$ligand,LR_i$receptor),,drop=FALSE]),
        check.names=FALSE
      )

      if(sum(exp_data[,1])==0 || sum(exp_data[,2])==0){

        b1 <- NA
        b2 <- NA
        b3 <- NA
        p2 <- NA
        p3 <- NA
        p_value <- NA

        x <- data.frame(LR_i,b1,b2,b3,p2,p3,p_value)
        TLR_model_res <- rbind(TLR_model_res,x)

      }else{

        Mutation_Status <- data.frame(
          t(matutation_0_1_1[LR_i$ligand,,drop=FALSE]),
          check.names=FALSE
        )

        ligandMutation <- exp_data[,1]*Mutation_Status[,1]

        data <- cbind(
          exp_data,
          Mutation_Status,
          ligandMutation,
          Tumor_Purity
        )

        colnames(data) <- c(
          "ligand",
          "receptor",
          "Mutation_Status",
          "ligandMutation",
          "Tumor_Purity"
        )

        model <- tryCatch(
          lm(
            receptor ~ ligand + Mutation_Status +
              ligandMutation + Tumor_Purity,
            data=data
          ),
          error=function(e) NULL
        )

        if(is.null(model)){

          b1 <- NA
          b2 <- NA
          b3 <- NA
          p2 <- NA
          p3 <- NA
          p_value <- NA

        }else{

          b1 <- model[["coefficients"]][["ligand"]]
          b2 <- model[["coefficients"]][["Mutation_Status"]]
          b3 <- model[["coefficients"]][["ligandMutation"]]

          if(is.na(b2) || is.na(b3)){

            p2 <- NA
            p3 <- NA
            p_value <- NA

          }else{

            p2 <- tryCatch(
              aod::wald.test(
                b=model[["coefficients"]],
                Sigma=vcov(model),
                Terms=3
              )[["result"]][["chi2"]][["P"]],
              error=function(e) NA_real_
            )

            p3 <- tryCatch(
              aod::wald.test(
                b=model[["coefficients"]],
                Sigma=vcov(model),
                Terms=4
              )[["result"]][["chi2"]][["P"]],
              error=function(e) NA_real_
            )

            if(is.finite(p2) && is.finite(p3)){
              p_values <- c(p2,p3)
              sum_log_p <- -2*sum(log(pmax(p_values,.Machine$double.xmin)))
              p_value <- pchisq(
                sum_log_p,
                df=2*length(p_values),
                lower.tail=FALSE
              )
            }else{
              p_value <- NA
            }
          }
        }

        x <- data.frame(LR_i,b1,b2,b3,p2,p3,p_value)
        TLR_model_res <- rbind(TLR_model_res,x)
      }

      if(i%%100==0){
        cat("TLR:",i,"/",nrow(TLR_unique_LR),"\n")
      }
    }
  }

  #####BH adjustment is performed separately for each direction at UNIQUE-LR level.
  #####Only successfully tested finite p-values enter each BH testing family.
  #####p_value retains Fisher(P2,P3); q is BH-adjusted Fisher p-value.
  oncoLR_BH <- function(p){
    adjusted <- rep(NA_real_,length(p))
    valid <- which(is.finite(p))
    if(length(valid)>0L){
      adjusted[valid] <- stats::p.adjust(p[valid],method="BH")
    }
    adjusted
  }

  TLR_model_res1 <- TLR_model_res[
    is.finite(TLR_model_res$p_value),,drop=FALSE
  ]

  if(nrow(TLR_model_res1)>0){
    TLR_model_res1$p_fisher <- TLR_model_res1$p_value
    TLR_model_res1$q2 <- oncoLR_BH(TLR_model_res1$p2)
    TLR_model_res1$q3 <- oncoLR_BH(TLR_model_res1$p3)
    TLR_model_res1$q_fisher <- oncoLR_BH(TLR_model_res1$p_fisher)
    TLR_model_res1$q <- TLR_model_res1$q_fisher

    TLR_model_res1$Condition <- ifelse(
      TLR_model_res1$q<q,
      "Sig",
      "No-sig"
    )

    TLR_model_res2 <- TLR_model_res1[
      TLR_model_res1$q<q,
      ,
      drop=FALSE
    ]
  }else{
    TLR_model_res1$p_fisher <- numeric(0)
    TLR_model_res1$q2 <- numeric(0)
    TLR_model_res1$q3 <- numeric(0)
    TLR_model_res1$q_fisher <- numeric(0)
    TLR_model_res1$q <- numeric(0)
    TLR_model_res1$Condition <- character(0)
    TLR_model_res2 <- TLR_model_res1
  }

  ########################################################################
  ##### LRT: Tumor is target; receptor mutation
  ########################################################################
  LRT_model_res <- data.frame()

  if(nrow(LRT_unique_LR)>0){
    for(i in seq_len(nrow(LRT_unique_LR))){

      LR_i <- LRT_unique_LR[i,,drop=FALSE]

      exp_data <- data.frame(
        t(exp12[c(LR_i$ligand,LR_i$receptor),,drop=FALSE]),
        check.names=FALSE
      )

      if(sum(exp_data[,1])==0 || sum(exp_data[,2])==0){

        b1 <- NA
        b2 <- NA
        b3 <- NA
        p2 <- NA
        p3 <- NA
        p_value <- NA

        x <- data.frame(LR_i,b1,b2,b3,p2,p3,p_value)
        LRT_model_res <- rbind(LRT_model_res,x)

      }else{

        Mutation_Status <- data.frame(
          t(matutation_0_1_1[LR_i$receptor,,drop=FALSE]),
          check.names=FALSE
        )

        receptorMutation <- exp_data[,2]*Mutation_Status[,1]

        data <- cbind(
          exp_data,
          Mutation_Status,
          receptorMutation,
          Tumor_Purity
        )

        colnames(data) <- c(
          "ligand",
          "receptor",
          "Mutation_Status",
          "receptorMutation",
          "Tumor_Purity"
        )

        model <- tryCatch(
          lm(
            ligand ~ receptor + Mutation_Status +
              receptorMutation + Tumor_Purity,
            data=data
          ),
          error=function(e) NULL
        )

        if(is.null(model)){

          b1 <- NA
          b2 <- NA
          b3 <- NA
          p2 <- NA
          p3 <- NA
          p_value <- NA

        }else{

          b1 <- model[["coefficients"]][["receptor"]]
          b2 <- model[["coefficients"]][["Mutation_Status"]]
          b3 <- model[["coefficients"]][["receptorMutation"]]

          if(is.na(b2) || is.na(b3)){

            p2 <- NA
            p3 <- NA
            p_value <- NA

          }else{

            p2 <- tryCatch(
              aod::wald.test(
                b=model[["coefficients"]],
                Sigma=vcov(model),
                Terms=3
              )[["result"]][["chi2"]][["P"]],
              error=function(e) NA_real_
            )

            p3 <- tryCatch(
              aod::wald.test(
                b=model[["coefficients"]],
                Sigma=vcov(model),
                Terms=4
              )[["result"]][["chi2"]][["P"]],
              error=function(e) NA_real_
            )

            if(is.finite(p2) && is.finite(p3)){
              p_values <- c(p2,p3)
              sum_log_p <- -2*sum(log(pmax(p_values,.Machine$double.xmin)))
              p_value <- pchisq(
                sum_log_p,
                df=2*length(p_values),
                lower.tail=FALSE
              )
            }else{
              p_value <- NA
            }
          }
        }

        x <- data.frame(LR_i,b1,b2,b3,p2,p3,p_value)
        LRT_model_res <- rbind(LRT_model_res,x)
      }

      if(i%%100==0){
        cat("LRT:",i,"/",nrow(LRT_unique_LR),"\n")
      }
    }
  }

  #####BH adjustment is performed independently from the TLR testing family.
  LRT_model_res1 <- LRT_model_res[
    is.finite(LRT_model_res$p_value),,drop=FALSE
  ]

  if(nrow(LRT_model_res1)>0){
    LRT_model_res1$p_fisher <- LRT_model_res1$p_value
    LRT_model_res1$q2 <- oncoLR_BH(LRT_model_res1$p2)
    LRT_model_res1$q3 <- oncoLR_BH(LRT_model_res1$p3)
    LRT_model_res1$q_fisher <- oncoLR_BH(LRT_model_res1$p_fisher)
    LRT_model_res1$q <- LRT_model_res1$q_fisher

    LRT_model_res1$Condition <- ifelse(
      LRT_model_res1$q<q,
      "Sig",
      "No-sig"
    )

    LRT_model_res2 <- LRT_model_res1[
      LRT_model_res1$q<q,
      ,
      drop=FALSE
    ]
  }else{
    LRT_model_res1$p_fisher <- numeric(0)
    LRT_model_res1$q2 <- numeric(0)
    LRT_model_res1$q3 <- numeric(0)
    LRT_model_res1$q_fisher <- numeric(0)
    LRT_model_res1$q <- numeric(0)
    LRT_model_res1$Condition <- character(0)
    LRT_model_res2 <- LRT_model_res1
  }

  ########################################################################
  ##### Step7: map significant UNIQUE LR models back to cell-resolved CCC
  ########################################################################
  print('Step7-------Map significant unique LR back to cell network')

  if(nrow(TLR_model_res2)>0){
    TLR_model_result <- merge(
      TCC_net_exp,
      TLR_model_res2,
      by=c("ligand","receptor")
    )

    TLR_model_result <- unique(TLR_model_result)

    TLR_first_cols <- c(
      "source",
      "ligand",
      "receptor",
      "target",
      "interaction_name"
    )

    TLR_model_result <- TLR_model_result[
      ,
      c(
        TLR_first_cols,
        setdiff(colnames(TLR_model_result),TLR_first_cols)
      ),
      drop=FALSE
    ]

  }else{

    TLR_model_result <- TCC_net_exp[0,,drop=FALSE]
  }

  if(nrow(LRT_model_res2)>0){
    LRT_model_result <- merge(
      CTC_net_exp,
      LRT_model_res2,
      by=c("ligand","receptor")
    )

    LRT_model_result <- unique(LRT_model_result)

    LRT_first_cols <- c(
      "source",
      "ligand",
      "receptor",
      "target",
      "interaction_name"
    )

    LRT_model_result <- LRT_model_result[
      ,
      c(
        LRT_first_cols,
        setdiff(colnames(LRT_model_result),LRT_first_cols)
      ),
      drop=FALSE
    ]

  }else{

    LRT_model_result <- CTC_net_exp[0,,drop=FALSE]
  }

  rownames(TLR_model_result) <- NULL
  rownames(LRT_model_result) <- NULL

  cat(
    "TLR significant unique LR:",
    nrow(TLR_model_res2),
    "; mapped cell events:",
    nrow(TLR_model_result),
    "\n"
  )

  cat(
    "LRT significant unique LR:",
    nrow(LRT_model_res2),
    "; mapped cell events:",
    nrow(LRT_model_result),
    "\n"
  )

  ############Return
  print('Step8-------Linear model is over')

  OncoLRI <- list(
    TLR_model_result,
    LRT_model_result
  )

  names(OncoLRI) <- c("TLR","LRT")

  #####Keep unique-LR model tables for QC without changing OncoLRI[[1]]/[[2]]
  attr(OncoLRI,"unique_LR_model_results") <- list(
    TLR=TLR_model_res1,
    LRT=LRT_model_res1
  )

  return(OncoLRI)
}

