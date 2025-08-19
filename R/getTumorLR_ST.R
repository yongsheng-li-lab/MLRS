#' @title Constructing a cell-cell interaction network involving tumor cells from the spatial transcriptome level.
#'
#' @param st_count A dataframe containing spatial transcriptome counts with each column representing a cell and each row representing a gene.
#' @param st_metadata A dataframe containing spatial transcriptomics coordinates with 4 columns, i.e., cell id, x, y indicating single cell spatial transcriptomics data and cell_types containing the ‘Tumor’ cell type.
#' @param ligand_receptor_DB Data.frame of ligands, receptors and interactions.
#' @param python_path path to python binary to use
#' @param outputDir Character, The output path of the currently running job where temporary and final results will be saved.
#'
#' @importFrom Giotto
#' @importFrom SpaTalk
#' @importFrom stMLnet
#' @importFrom dplyr
#' @importFrom tidyverse
#'
#' @return
#' @export
#'
#' @examples
getTumorLR_ST <- function(st_count,st_metadata,ligand_receptor_DB,python_path,outputDir){
  #########Run Giotto
  print('Giotto is running')
  myinstructions = createGiottoInstructions(save_dir = NULL,
                                            save_plot = FALSE,
                                            show_plot = F,
                                            python_path = python_path)
  spot_loc <- st_metadata[,c("spot","imagerow","imagecol")]
  colnames(spot_loc) <- c('cell_ID','sdimx','sdimy')
  visium <- createGiottoObject(raw_exprs = st_count,
                               spatial_locs = spot_loc,
                               instructions = myinstructions)
  visium <- normalizeGiotto(gobject = visium, scalefactor = 10000, verbose = T)
  visium <- addStatistics(gobject = visium)
  visium <- calculateHVG(gobject = visium)
  visium <- runPCA(gobject = visium)
  visium <- runUMAP(visium, dimensions_to_use = 1:30)
  visium <- runtSNE(visium, dimensions_to_use = 1:30,check_duplicates = FALSE)
  visium <- createNearestNetwork(gobject = visium, dimensions_to_use = 1:30, k = 15)
  visium <- doLeidenCluster(gobject = visium, resolution = 0.4, n_iterations = 1000)
  visium@cell_metadata[["cell_ID"]] <- st_metadata$spot
  visium@cell_metadata$cell_type <- st_metadata$predict_celltype
  Giotto_LR_data <- ligand_receptor_DB
  Giotto_LR_data$ligand_det <- ifelse(Giotto_LR_data$ligand %in% visium@gene_ID, T, F)
  Giotto_LR_data$receptor_det <- ifelse(Giotto_LR_data$receptor %in% visium@gene_ID, T, F)
  Giotto_LR_data_det = Giotto_LR_data[which(Giotto_LR_data$ligand_det == T & Giotto_LR_data$receptor_det == T),]
  select_ligands = Giotto_LR_data_det$ligand
  select_receptors = Giotto_LR_data_det$receptor
  expr_only_scores = exprCellCellcom(gobject = visium,
                                     cluster_column = 'cell_type',
                                     random_iter = 500,
                                     gene_set_1 = select_ligands,
                                     gene_set_2 = select_receptors)
  visium = createSpatialNetwork(gobject = visium, minimum_k = 2, maximum_distance_delaunay = 400)
  spatial_all_scores = spatCellCellcom(visium,
                                       spatial_network_name = 'Delaunay_network',
                                       cluster_column = 'cell_type',
                                       random_iter = 500,
                                       gene_set_1 = select_ligands,
                                       gene_set_2 = select_receptors,
                                       adjust_method = 'fdr',
                                       do_parallel = T,
                                       cores = 4,
                                       verbose = 'none')
  comb_comm = combCCcom(spatialCC = spatial_all_scores,
                        exprCC = expr_only_scores)
  Gitto_res <- comb_comm[-which(comb_comm$lig_expr==0),]
  Gitto_res <- Gitto_res[-which(Gitto_res$rec_expr==0),]
  Gitto_res <- Gitto_res[-which(Gitto_res$lig_expr_spat==0),]
  Gitto_res <- Gitto_res[-which(Gitto_res$rec_expr_spat==0),]
  Gitto_net <- Gitto_res[,c("lig_cell_type","ligand","receptor","rec_cell_type")]
  Gitto_net$interaction_name <- paste0(Gitto_net$ligand,"_",Gitto_net$receptor)
  colnames(Gitto_net) <- c("source","ligand","receptor","target","interaction_name")
  print('Giotto over')
  #########Run SpaTalk
  print('SpaTalk is running')
  spot_loc <- st_metadata[,c("spot","imagerow","imagecol")]
  colnames(spot_loc) <- c('cell','x','y')
  SpaTalk_data <- createSpaTalk(st_count,
                                spot_loc,
                                'Human',
                                if_st_is_sc=T,
                                celltype=st_metadata$predict_celltype)
  SpaTalk_Giotto_LR_data <- ligand_receptor_DB
  SpaTalk_Giotto_LR_data$species='Human'
  SpaTalk_obj <- find_lr_path(object = SpaTalk_data,
                              lrpairs = SpaTalk_Giotto_LR_data,
                              pathways = pathways)
  SpaTalk_obj <- dec_cci_all(SpaTalk_obj)
  SpaTalk_net <- SpaTalk_obj@lrpair[,c("celltype_sender","ligand","receptor","celltype_receiver")]

  SpaTalk_net$interaction_name <- paste0(SpaTalk_net$ligand,"_",SpaTalk_net$receptor)
  colnames(SpaTalk_net) <- c("source","ligand","receptor","target","interaction_name")
  print('SpaTalk over')
  #########Run stMLnet_net
  print('stMLnet_net is running')
  spot_loc <- st_metadata[,c("spot","imagerow","imagecol")]
  colnames(spot_loc) <- c('cell','x','y')
  spot_cell_type <- st_metadata[,c("spot","predict_celltype")]
  colnames(spot_cell_type) <- c('Barcode','Cluster')
  outputDir <- outputDir
  data('stMLnet_ex_databases',package = 'MLRS')
  stMLnet_net_LR_data <- ligand_receptor_DB
  colnames(stMLnet_net_LR_data) <- c('source','target')
  stMLnet_ex_databases$LigRec.DB <- stMLnet_net_LR_data
  ex_inputs <- MLRS::Select_Lig_Rec_TGs(ExprMat = st_count,
                                        AnnoMat = spot_cell_type,
                                        LocaMat = spot_loc,
                                        Databases = stMLnet_ex_databases,
                                        python_path = python_path,
                                        min.pct = 0.05,
                                        expr.ct = 0.1,
                                        pct.ct = 0.05)
  restMLnet <- runMLnet(ExprMat = ex_inputs$exprMat, AnnoMat = ex_inputs$annoMat,
                        LigClus = NULL, RecClus = NULL, Normalize = F,
                        OutputDir = outputDir, Databases = stMLnet_ex_databases,
                        TGList=ex_inputs$tgs_of_inter, LigList=ex_inputs$ligs_of_inter, RecList=ex_inputs$recs_of_inter)
  file_list <- list.files(pattern = "cellpair_LRI_pval.rds",recursive = T)
  data_stMLnet_net <- data.frame()
  for(i in 1:length(file_list)){
    data <- readRDS(file_list[i])
    if(nrow(data)>0){
      data$cell_cell <- strsplit(file_list[i],"/")[[1]][2]
      data_stMLnet_net <- rbind(data_stMLnet_net,data)
      print(i)
    }
  }
  stMLnet_celltype <- matrix(unlist(strsplit(data_stMLnet_net$cell_cell,"_")),nrow = length(data_stMLnet_net$cell_cell),byrow = T)
  colnames(stMLnet_celltype) <- c('source_cell','target_cell')
  stMLnet_net <- data.frame(data_stMLnet_net[,c(1,2)],stMLnet_celltype)
  stMLnet_net <- stMLnet_net[,c("source_cell","source","target","target_cell")]
  stMLnet_net$interaction_name <- paste0(stMLnet_net$source,"_",stMLnet_net$target)
  colnames(stMLnet_net) <- c("source","ligand","receptor","target","interaction_name")
  print('stMLnet over')
  #######Return
  Gitto_net$Gitto <- 1
  stMLnet_net$stMLnet <- 1
  SpaTalk_net$SpaTalk <- 1
  Gitto_stMLnet_net <- full_join(Gitto_net,stMLnet_net,by=c("source","ligand","receptor","target","interaction_name"),relationship = "many-to-many")
  Gitto_stMLnet_SpaTalk_net <- full_join(Gitto_stMLnet_net,SpaTalk_net,by=c("source","ligand","receptor","target","interaction_name"),relationship = "many-to-many")
  Gitto_stMLnet_SpaTalk_net[is.na(Gitto_stMLnet_SpaTalk_net)] <- 0
  STCCC_net <- Gitto_stMLnet_SpaTalk_net
  STCCC_net$Freq <- rowSums(STCCC_net[,6:ncol(STCCC_net)])
  return(STCCC_net)
}
