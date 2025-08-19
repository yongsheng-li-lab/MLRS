#' @title Organizing CCC results for multiple idle slices based on the dplyr package.
#'
#' @param net1 CCC results for slice 1 based on inferences from the getTumorLR_ST() function.
#' @param net2 CCC results for slice 2 based on inferences from the getTumorLR_ST() function.
#' @param is.union Default TRUE. Take the union set of the CCC results inferred for the slices.
#' @param is.intersect Default FALSE，if TRUE，Take the intersect set of the CCC results inferred for the slices.
#'
#' @importFrom dplyr full_join inner_join
#'
#' @return
#' @export
#'
#' @examples
TumorLR_multiSlide <- function(net1,net2,is.union=T,is.intersect=F){
  if(is.union==T&is.intersect==F){
    slide_STCCC_net <- full_join(net1,net2,by=colnames(net1),relationship = "many-to-many")
    return(slide_STCCC_net)
  }
  if(is.union==F&is.intersect==T){
    slide_STCCC_net <- inner_join(net1,net2,by=colnames(net1),relationship = "many-to-many")
    return(slide_STCCC_net)
  }
  if(is.union==T&is.intersect==T){
    slide_STCCC_net_union <- full_join(net1,net2,by=colnames(net1),relationship = "many-to-many")
    slide_STCCC_net2_intersect <- inner_join(net1,net2,by=colnames(net1),relationship = "many-to-many")
    return(list(slide_STCCC_net_union,slide_STCCC_net2_intersect))
  }
}
