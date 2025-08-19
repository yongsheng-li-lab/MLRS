#' @title Calculate mutation frequency
#'
#' @param mutation_obj a dataframe generated from the TCGA-MC3 project data, the user-defined dataset need to include at least Hugo_Symbol, Tumor_Sample_Barcode, Variant_Classification columns.
#' @importFrom maftools plotmafSummary plotTiTv
#' @importFrom utils read.delim
#'
#' @return mutation_obj_frequency
#'
#' @export
#'
#' @examples
mutfreq <- function(mutation_obj){
  #############################calculate cancer mutation##############################
  cancer_gene <- as.character(unique(mutation_obj$Hugo_Symbol))
  cancer_sample <- as.character(unique(mutation_obj$Tumor_Sample_Barcode))
  cancer_matutation <- as.data.frame(matrix("",length(cancer_gene),length(cancer_sample),
                                               dimnames = list(cancer_gene,cancer_sample)))
  cancer_matutation_0_1 <- as.data.frame(matrix(0,length(cancer_gene),length(cancer_sample),
                                                   dimnames = list(cancer_gene,cancer_sample)))
  for (i in 1:nrow(mutation_obj)){
    cancer_matutation[as.character(mutation_obj[i,"Hugo_Symbol"]),as.character(mutation_obj[i,"Tumor_Sample_Barcode"])] <- as.character(mutation_obj[i,"Variant_Classification"])
    print(i)
  }#Show what types of mutations occurred in each gene in each sample
  for (i in 1:nrow(mutation_obj)){
    cancer_matutation_0_1[as.character(mutation_obj[i,"Hugo_Symbol"]),as.character(mutation_obj[i,"Tumor_Sample_Barcode"])] <- 1
    print(i)
  }#Show which mutation types occurred in each gene in each sample, here the mutation type is mutated directly to 1 and not 0
  #All samples mutations summarized/ranked
  cancer_gene_count <- data.frame(cancer_gene=rownames(cancer_matutation_0_1),
                                     count=as.numeric(apply(cancer_matutation_0_1,1,sum))) %>%arrange(desc(count))
  mutation_obj_frequency <- cancer_gene_count
  mutation_obj_frequency$cancer_sample_sum <- c(rep(length(cancer_sample),nrow(mutation_obj_frequency)))
  mutation_obj_frequency$frequency <- round(mutation_obj_frequency$count/length(cancer_sample),4)
  colnames(mutation_obj_frequency) <- c('gene','count','sample_sum','frequency')
  return(mutation_obj_frequency)
}
