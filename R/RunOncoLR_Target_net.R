#' @title Constructs cell-to-cell post-preference Ligand-Receptor-TarGene networks based on scRNA-Seq expression matrices and cell_type tables
#'
#' @param sc_data A Seurat object
#' @param cell_type A character vector representing ligand/recepor cell types.
#' @param oncoLR A dataframe, including ligand, receptor, interaction_name.
#' @param ligand.source.cell Ligand source cells. Default 'Tumor cell'
#' @param receptor.source.cell Receptor source cells. Default 'Other cell'
#' @param ligand_target_matrix Indicates the regulation score between ligand-target (column is ligand). The default setting is the file “ligand_target_matrix_nsga2r_final.rds” in the /databases/ folder.
#' @param weighted_networks Weighted network representing interactions and their weights in the ligand signaling and gene regulatory networks. The default setting is the file “weighted_networks_nsga2r_final.rds” in the /databases/ folder.
#'
#' @importFrom nichenetr
#'
#' @return
#' @export
#'
#' @examples
RunOncoLR_Target_net <- function(sc_data,cell_type,oncoLR,ligand.source.cell='Tumor cell',receptor.source.cell='Other cell',ligand_target_matrix,weighted_networks){
  if(ligand.source.cell=='Tumor cell'){
    oncoLR_target_res <- data.frame()
    for(i in cell_type){
      TLR_network <- oncoLR[which(oncoLR$target==i),c("ligand","receptor")]
      colnames(TLR_network) <- c('from','to')
      TLR_network_ligand <- data.frame(intersect(TLR_network$from,colnames(ligand_target_matrix)))
      colnames(TLR_network_ligand) <- c('from')
      TLR_network1 <- merge(TLR_network,TLR_network_ligand,by='from')
      ligand_target_matrix1 <- ligand_target_matrix[,TLR_network_ligand$from]
      nichenet_output1 <- nichenet_seuratobj_cluster_de(
                          seurat_obj = sc_data,
                          receiver_affected = i,
                          receiver_reference = "Tumor",
                          sender = c("Tumor"),
                          ligand_target_matrix = ligand_target_matrix1,
                          lr_network = TLR_network1,
                          weighted_networks = weighted_networks)
      LR <- nichenet_output1[["ligand_receptor_df"]]
      LTarget <- nichenet_output1[["ligand_target_df"]]
      LR_target <- merge(LR,LTarget,by='ligand')
      source <- data.frame(rep("Tumor",nrow(LR_target)))
      target <- data.frame(rep(i,nrow(LR_target)))
      x <- cbind(source,LR_target,target)
      oncoLR_target_res <- rbind(oncoLR_target_res,x)
      print(i)
    }
    colnames(oncoLR_target_res) <- c('source','ligand','receptor','LR_weight','Target_gene','LTarget_weight','target')
    return(oncoLR_target_res)
  }
  if(receptor.source.cell=='Tumor cell'){
    oncoLR_target_res <- data.frame()
    for(i in cell_type){
      LRT_network <- oncoLR[which(oncoLR$source==i),c("ligand","receptor")]
      colnames(LRT_network) <- c('from','to')
      LRT_network_ligand <- data.frame(intersect(LRT_network$from,colnames(ligand_target_matrix)))
      colnames(LRT_network_ligand) <- c('from')
      LRT_network1 <- merge(LRT_network,LRT_network_ligand,by='from')
      ligand_target_matrix2 <- ligand_target_matrix[,LRT_network_ligand$from]
      nichenet_output1 <- nichenet_seuratobj_cluster_de(
        seurat_obj = sc_data,
        receiver_affected = "Tumor",
        receiver_reference = i,
        sender =i,
        ligand_target_matrix = ligand_target_matrix2,
        lr_network = LRT_network1,
        weighted_networks = weighted_networks)
      LR <- nichenet_output1[["ligand_receptor_df"]]
      LTarget <- nichenet_output1[["ligand_target_df"]]
      LR_target <- merge(LR,LTarget,by='ligand')
      source <- data.frame(rep(i,nrow(LR_target)))
      target <- data.frame(rep("Tumor",nrow(LR_target)))
      x <- cbind(source,LR_target,target)
      oncoLR_target_res <- rbind(oncoLR_target_res,x)
      print(i)
    }
    colnames(oncoLR_target_res) <- c('source','ligand','receptor','LR_weight','Target_gene','LTarget_weight','target')
    return(oncoLR_target_res)
  }
}
