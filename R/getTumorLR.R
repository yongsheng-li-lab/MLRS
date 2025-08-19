#' @title Organizing CCC results for single-cell data and multiple idle slices based on the dplyr package.
#'
#' @param SC_net CCC results inferred from single-cell data.
#' @param ST_net CCC results inferred from spatial transcriptome data.
#' @param is.union Default TRUE. Take the union set of the CCC results inferred for the slices.
#' @param is.intersect Default FALSE，if TRUE，Take the intersect set of the CCC results inferred for the slices.
#'
#' @importFrom dplyr full_join inner_join
#'
#' @return
#' @export
#'
#' @examples
getTumorLR <- function(SC_net,ST_net,is.union=T,is.intersect=F){
  ######Read data
  SC_net <- SC_net[,-ncol(SC_net)]
  ST_net <- ST_net[,-ncol(ST_net)]
  ######Run
  if(is.union==T&is.intersect==F){
    SC_ST_CCC_net <- full_join(SC_net,ST_net,by=colnames(SC_net)[1:5],relationship = "many-to-many")
    SC_ST_CCC_net[is.na(SC_ST_CCC_net)] <- 0
    SC_ST_CCC_net$Freq <- rowSums(SC_ST_CCC_net[,6:ncol(SC_ST_CCC_net)])
    return(SC_ST_CCC_net)
  }
  if(is.union==F&is.intersect==T){
    SC_ST_CCC_net <- inner_join(SC_net,ST_net,by=colnames(SC_net)[1:5],relationship = "many-to-many")
    SC_ST_CCC_net[is.na(SC_ST_CCC_net)] <- 0
    SC_ST_CCC_net$Freq <- rowSums(SC_ST_CCC_net[,6:ncol(SC_ST_CCC_net)])
    return(SC_ST_CCC_net)
  }
  if(is.union==T&is.intersect==T){
    SC_ST_CCC_net_union <- full_join(SC_net,ST_net,by=colnames(SC_net)[1:5],relationship = "many-to-many")
    SC_ST_CCC_net_union[is.na(SC_ST_CCC_net_union)] <- 0
    SC_ST_CCC_net_union$Freq <- rowSums(SC_ST_CCC_net_union[,6:ncol(SC_ST_CCC_net_union)])
    SC_ST_CCC_net_intersect <- inner_join(SC_net,ST_net,by=colnames(SC_net)[1:5],relationship = "many-to-many")
    SC_ST_CCC_net_intersect[is.na(SC_ST_CCC_net_intersect)] <- 0
    SC_ST_CCC_net_intersect$Freq <- rowSums(SC_ST_CCC_net_intersect[,6:ncol(SC_ST_CCC_net_intersect)])
    return(list(SC_ST_CCC_net_union,SC_ST_CCC_net_intersect))
  }
}
