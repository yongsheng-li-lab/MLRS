#' @title Identify the methods in which the Ligand-Receptor is activated.
#'
#' @param CCC_net Cell-cell communication network.
#' @param LR_pairs A String of ligand-receptor connected by “_”.
#'
#' @importFrom dplyr
#'
#' @return
#' @export
#'
#' @examples
getLR_activate <- function(CCC_net,LR_pairs){
  if('interaction_name' %in% colnames(CCC_net)){
    LR_pairs <- data.frame(LR_pairs)
    colnames(LR_pairs) <- c('interaction_name')
    LR_activate <- inner_join(CCC_net,LR_pairs,by='interaction_name')
  }else{
    CCC_net$interaction_name <- paste0(CCC_net$ligand,"_",CCC_net$receptor)
    LR_pairs <- data.frame(LR_pairs)
    colnames(LR_pairs) <- c('interaction_name')
    LR_activate <- inner_join(CCC_net,LR_pairs,by='interaction_name')
  }

}
