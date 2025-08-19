#' @title Constructing a cell-cell interaction network involving tumor cells from the single-cell level.
#'
#' @param sc_data Seurat object to create ligand receptor matrices from.
#' @param sc_metadata Grouping variable in the Seurat object metadata used to define the groups of cells. Defaults to two columns, cell_id and cell_types containing the ‘Tumor’ cell type.
#' @param ligand_receptor_DB Data.frame of ligands, receptors and interactions，The default is the TICC and Omnipath datasets. If NULL, database of ligand/receptors interactions that comes with the method is used.
#' @param sc_count Pre-processed (i.e. no raw counts) gene expression matrix with genes on rows and cells on columns; rownames(gene_expr) must contain the genes names.
#'
#' @importFrom Seurat
#' @importFrom SeuratData
#' @importFrom CellChat
#' @importFrom celltalker
#' @importFrom SIMLR
#' @importFrom SingleCellSignalR
#' @importFrom iTALK
#' @importFrom icellnet
#' @importFrom stringr
#'
#' @return
#' @export
#'
#' @examples
getTumorLR_SC <- function(sc_data,sc_count,sc_metadata,ligand_receptor_DB=NULL){
  #####Reading sc data and its metadata
  Idents(sc_data) <- sc_metadata$cell_type_final
  sc_data@meta.data$cell_type_final= Idents(sc_data)
  if(is.null(ligand_receptor_DB)!=T){
    #####Reading ligand_receptor_DB
    LR_data <- ligand_receptor_DB
    LR_data <- LR_data[!duplicated(LR_data),]
    ligand <- unique(LR_data$ligand)
    receptor <- unique(LR_data$receptor)
    db <- LR_data
    colnames(db) <- c("ligand","receptor")
    db$interaction_name <- paste0(db$ligand,"_",db$receptor)
    #####Run cellchat
    print('CellChat is running')
    cellchat <- createCellChat(object = sc_data, meta = sc_data@meta.data, group.by = "cell_type_final" )
    gene_info <- data.frame(unique(c(db$ligand,db$receptor)))
    colnames(gene_info) <- c("Symbol")
    gene_info$EntryName <- paste0(gene_info$Symbol,"_","HUMAN")
    cellchatDB <- updateCellChatDB(db = db, gene_info = gene_info)
    cellchat@DB <- cellchatDB
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat_net <- subsetCommunication(cellchat)
    cellchat_net <- cellchat_net[,c("source","ligand","receptor","target","interaction_name")]
    print('CellChat is over')
    #####Run scSeqComm
    print('scSeqComm is running')
    gene_expr_matrix <- as.matrix(sc_count)
    cell_metadata <- data.frame(colnames(as.matrix(sc_count)),sc_metadata$cell_type_final)
    colnames(cell_metadata) <- c("Cell_ID","Cluster_ID")
    scSeqComm <- scSeqComm_analyze(gene_expr_matrix,
                                        cell_metadata = cell_metadata,
                                        inter_signaling = TRUE,
                                        LR_pairs_DB = db,
                                        intra_signaling = FALSE,
                                        N_cores = 1)
    rm(gene_expr_matrix)
    scSeqComm_result <- scSeqComm[["comm_results"]]
    scSeqComm_result <- scSeqComm_result[-which(scSeqComm_result$mean.expression==0),]
    scSeqComm_net <- scSeqComm_result[,c("cluster_L","ligand","receptor","cluster_R")]
    scSeqComm_net$interaction_name <- paste0(scSeqComm_net$ligand,"_",scSeqComm_net$receptor)
    colnames(scSeqComm_net) <- c("source","ligand","receptor","target","interaction_name")
    print('scSeqComm is over')
    #####Run cellTalker
    print('cellTalker is running')
    cellTalker <- celltalk(input_object=sc_data,
                            metadata_grouping="cell_type_final",
                            ligand_receptor_pairs=db,
                            number_cells_required=100,
                            min_expression=1000,
                            max_expression=20000,
                            scramble_times=10)
    cellTalker_net <- cellTalker[,c("cell_type1","cell_type2","interaction")]
    cellTalker_net_LR <- matrix(unlist(strsplit(cellTalker_net$interaction,"_")),nrow = length(cellTalker_net$interaction),byrow = T)
    colnames(cellTalker_net_LR) <- c("ligand",'receptor')
    cellTalker_net_LR <- as.data.frame(cellTalker_net_LR)
    cellTalker_net <- cbind(cellTalker_net,cellTalker_net_LR)
    cellTalker_net <- cellTalker_net[,c("cell_type1","ligand",'receptor',"cell_type2","interaction")]
    colnames(cellTalker_net) <- c("source","ligand","receptor","target","interaction_name")
    print('cellTalker is over')
    #####Run Connectome
    print('Connectome is running')
    Connectome <- CreateConnectome(
      sc_data,
      LR.database = "custom",
      species = 'human',
      custom.list =db,
      min.cells.per.ident=75
    )
    Connectome <- Connectome[-which(Connectome$recept.expression==0),]
    Connectome <- Connectome[-which(Connectome$ligand.expression==0),]
    Connectome_net <- Connectome[,c("source","ligand","receptor","target")]
    Connectome_net$interaction_name <- paste0(Connectome_net$ligand,"_",Connectome_net$receptor)
    print('Connectome is over')
    #####Run iTALK
    print('iTALK is running')
    gene_expr <- as.data.frame(t(as.matrix(sc_data@assays$RNA@counts)))
    gene_expr$cell_type <- sc_data$cell_type_final
    gene_expr$compare_group <- sc_data$cell_type_final
    highly_gene_exprrs_genes<-rawParse(gene_expr,top_genes=50,stats='mean')
    iTALK_db <- data.frame(db$interaction_name,db$ligand,db$ligand,db$receptor,db$receptor,rep('other',nrow(db)))
    colnames(iTALK_db) <- colnames(iTALK:::database)
    res_cat<-FindLR(highly_gene_exprrs_genes,datatype='mean count',comm_type="other",database=iTALK_db)
    iTALK<-res_cat[order(res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs,decreasing=T),]
    iTALK<-iTALK[order(iTALK$cell_from_mean_exprs*iTALK$cell_to_mean_exprs,decreasing=T),]
    rm(gene_expr)
    iTALK <- iTALK[iTALK[["cell_from_mean_exprs"]] != 0, ]
    iTALK <- iTALK[iTALK[["cell_to_mean_exprs"]] != 0, ]
    iTALK_net <- iTALK[,c("cell_from","ligand","receptor","cell_to")]
    colnames(iTALK_net) <- c("source","ligand","receptor","target")
    iTALK_net$interaction_name <- paste0(iTALK_net$ligand,"_",iTALK_net$receptor)
    print('iTALK is over')
    #####Run icellnet
    print('icellnet is running')
    ICELLNETdb <- data.frame(cbind(db$ligand,db$ligand,db$ligand,db$ligand,db$receptor,db$receptor,db$receptor,db$receptor,db$receptor))
    colnames(ICELLNETdb) <- c('Ligand 1','Ligand 2','Ligand 3','Ligand 4','Receptor 1','Receptor 2','Receptor 3','Receptor 4','Receptor 5')
    ICELLNETdb[,c('Ligand 2','Ligand 3','Ligand 4','Receptor 2','Receptor 3','Receptor 4','Receptor 5')] <- NA
    ICELLNETdb$Family <- NA
    ICELLNETdb$Subfamily <- NA
    icellnet_average_clean= sc.data.cleaning(object = sc_data, db = ICELLNETdb, filter.perc = 0, save_file = F, path="", force.file = F)
    icellnet_data=as.data.frame(gene.scaling(as.data.frame(icellnet_average_clean), n=1, db=ICELLNETdb))
    icellnet_PC.target=data.frame("Class"=unique(sc_metadata$cell_type_final), "ID"= unique(sc_metadata$cell_type_final), "Cell_type"=unique(sc_metadata$cell_type_final))
    rownames(icellnet_PC.target)=icellnet_PC.target$Cell_type
    icellnet <- NULL
    for(i in icellnet_PC.target$Cell_type){
      score.computation= icellnet.score(direction="out", PC.data=icellnet_data,
                                        CC.data= as.data.frame(icellnet_data[,i], row.names = rownames(icellnet_data)),
                                        PC.target = icellnet_PC.target, PC=icellnet_PC.target$Cell_type, CC.type = "RNAseq",
                                        PC.type = "RNAseq",db=ICELLNETdb)
      LR_score=data.frame(score.computation[[2]])
      colnames(LR_score) <- icellnet_PC.target$Cell_type
      LR_score$cell_type <- i
      icellnet <- rbind(icellnet,LR_score)
      score.computation <- c()
      print(i)
    }
    icellnet1 <- icellnet
    icellnet1$ligand <- ICELLNETdb$`Ligand 1`
    icellnet1$receptor <- ICELLNETdb$`Receptor 1`
    colnames(icellnet1) <- c(icellnet_PC.target$Cell_type,'source_cell_type',"ligand",'receptor')
    icellnet_net <- NULL
    for(i in icellnet_PC.target$Cell_type){
      x <- data.frame(icellnet1[,c("source_cell_type","ligand","receptor")],icellnet1[,i])
      x$target_cell_type <- i
      icellnet_net <- rbind(icellnet_net,x)
      x <- c()
      print(i)
    }
    colnames(icellnet_net) <- c("source_cell_type","ligand","receptor",'LR_score','target_cell_type')
    icellnet_net <- icellnet_net[,c("source_cell_type","ligand","receptor",'target_cell_type','LR_score')]
    icellnet_net <- na.omit(icellnet_net)
    icellnet_net <- icellnet_net[-which(icellnet_net$LR_score==0),]
    rownames(icellnet_net) <- 1:nrow(icellnet_net)
    icellnet_net$interaction_name <- paste0(icellnet_net$ligand,"_",icellnet_net$receptor)
    ICELLNET_net <- icellnet_net[,c("source_cell_type","ligand","receptor","target_cell_type","interaction_name")]
    colnames(ICELLNET_net) <- c("source","ligand","receptor","target","interaction_name")
    print('icellnet is over')
    #####Run SingleCellSignalR
    print('SingleCellSignalR is running')
    gene_expr_dataframe <- as.data.frame(sc_count)
    cluster <- data.frame(sc_data$cell_type_final)
    colnames(cluster) <- c('cell_type_final')
    cell_type_cluster <- data.frame(unique(sc_data$cell_type_final))
    colnames(cell_type_cluster) <- c('cell_type_final')
    cell_type_cluster$cluster <- 1:nrow(cell_type_cluster)
    cluster1 <- inner_join(cluster,cell_type_cluster,by = join_by(cell_type_final))
    cluster1$cluster <- as.numeric(cluster1$cluster)
    SingleCellSignalR_cell_signaling <- function(data, genes,
                                                 cluster,int.type=c("paracrine","autocrine"),
                                                 c.names=NULL,s.score=0.5,logFC=log2(1.5),
                                                 species=c("homo sapiens","mus musculus"),
                                                 tol=0,write=TRUE,verbose=TRUE){
      if (dir.exists("cell-signaling")==FALSE & write==TRUE){
        dir.create("cell-signaling")
      }
      if (is.null(c.names)==TRUE){
        c.names <- paste("cluster",seq_len(max(cluster)))
      }
      if (min(cluster)!=1){
        cluster <- cluster + 1 - min(cluster)
      }
      if (length(c.names)!=max(cluster) | sum(duplicated(c.names))>0 |
          grepl("/",paste(c.names,collapse =""))){
        stop("The length of c.names must be equal to the number of clusters
        and must contain no duplicates. The cluster names must not include
        special characters")
      }
      int.type <- match.arg(int.type)
      species <- match.arg(species)

      rownames(data) <- genes
      z <- seq_len(max(cluster))
      lig <- unique(SingleCellSignalR_LRDB$ligand)
      rec <- unique(SingleCellSignalR_LRDB$receptor)
      data <- data.frame(data)
      data <- data[rowSums(data)>0,]
      med <- sum(data)/(nrow(data)*ncol(data))

      if (species=='mus musculus'){
        Hs2mm <- mm2Hs[,1]
        mm2Hs <- mm2Hs[,2]
        names(mm2Hs) <- as.character(Hs2mm)
        names(Hs2mm) <- as.character(mm2Hs)
        m.names <- mm2Hs[rownames(data)]
        data <- subset(data,(!is.na(m.names)))
        m.names <- m.names[!is.na(m.names)]
        rownames(data) <- as.character(m.names)
      }
      if (int.type=="autocrine"){
        if (verbose==TRUE){
          cat("Autocrine signaling: ",fill=TRUE)
        }
        auto <- list()
        k=0
        int=NULL
        n.int=NULL
        if (verbose==TRUE){
          cat("Checking for cell/cell signaling:",fill=TRUE)
        }
        for (i in z){
          if (sum(cluster==i)>1){
            tmp <- data[,cluster==i]
            tmp <- tmp[rowSums(tmp)>0,]
            if (sum(is.element(lig, rownames(tmp)))>0){
              lig.tmp <- rownames(tmp)[is.element(rownames(tmp),lig)]
              #lig.tmp <- lig.tmp[!is.element(lig.tmp,gene.list[[i]])]
            } else {lig.tmp=NULL}

            final.tmp <- SingleCellSignalR_LRDB[is.element(SingleCellSignalR_LRDB$ligand,lig.tmp),seq_len(2)]
            final.tmp <- data.frame(final.tmp,as.character(
              rep("autocrine|paracrine",sum(is.element(SingleCellSignalR_LRDB$ligand,lig.tmp)))))
            m.lig <- rowSums(tmp[unique(final.tmp[,1]),])/sum(cluster==i)
            names(m.lig) <- unique(final.tmp[,1])

            if (sum(is.element(rec, rownames(tmp)))>0){
              rec.tmp <- rownames(tmp)[is.element(rownames(tmp[apply(
                tmp,1,function(x) sum(x>0))>tol*ncol(tmp),]),rec)]
            } else {rec.tmp=NULL}

            for (j in z){
              if (sum(cluster==i)>1){
                temp <- data[,cluster==j]
                temp <- temp[rowSums(temp)>0,]
                if (sum(is.element(rec, rownames(temp)))>0){
                  rec.temp <- rownames(temp)[is.element(rownames(temp),rec)]
                  #rec.temp <- rec.temp[!is.element(rec.temp,gene.list[[j]])]
                } else {rec.temp=NULL}
                rec.temp <- rec.temp[is.element(rec.temp,rec.tmp)]
                m.rec <- rowSums(data.frame(temp[rec.temp,]))/sum(cluster==j)
                names(m.rec) <- rec.temp

                final <- final.tmp[is.element(final.tmp$receptor,rec.temp),]
                final <- cbind(final,SingleCellSignalR_LRscore(m.lig[final$ligand],m.rec[final$receptor],
                                                               med))

                colnames(final) <- c(c.names[i],c.names[j],"interaction type","SingleCellSignalR_LRscore")

                if (i==j){
                  final$`interaction type`="autocrine"
                }

                final <- final[final[,4]>s.score,]
                final <- final[order(final[,4],decreasing=TRUE),]

                if (species=="mus musculus"){
                  final[,1] <- Hs2mm[as.character(final[,1])]
                  final[,2] <- Hs2mm[as.character(final[,2])]
                }

                if (nrow(final)>0){
                  k <- k+1
                  auto[[k]] <- final
                  if (verbose==TRUE){
                    cat(paste(nrow(final),"interactions from",c.names[i],
                              "to",c.names[j]),fill=TRUE)
                  }
                  int <- c(int,paste(i,"-",j,sep=""))
                  n.int <- c(n.int,paste(c.names[i],"-",c.names[j],sep=""))
                  gr <- graph_from_data_frame(final,directed=FALSE)
                  if (write==TRUE){
                    fwrite(data.frame(final),paste("./cell-signaling/LR_interactions_",
                                                   c.names[i],"-",c.names[j],"-",
                                                   int.type,".txt",sep=""),sep="\t")
                  }
                }
              }
            }
          }
        }
        if (k!=0){
          names(auto) <- n.int
        }
      }


      ## Paracrine -------------------
      if (int.type=="paracrine"){
        gene.list <- vector("list",max(cluster))
        for (i in z){
          if (file.exists(paste("./cluster-analysis/table_dge_",c.names[i],
                                ".txt",sep=""))==TRUE){
            resu <- fread(paste("./cluster-analysis/table_dge_",c.names[i],
                                ".txt",sep=""),data.table=FALSE)
            gene.list[[i]] <- resu$genes[resu$logFC>logFC]
            if (species == "mus musculus"){
              gene.list[[i]] <- mm2Hs[gene.list[[i]]]
            }
            gene.list[[i]] <- gene.list[[i]][!is.na(gene.list[[i]])]
          } else {
            gene.list[[i]] <- "none"
            cat(paste("No such file as table_dge_",c.names[i],
                      ".txt in the cluster-analysis folder", sep=""),fill=TRUE)
          }
        }

        if (verbose==TRUE){
          cat("Paracrine signaling: ",fill=TRUE)
        }
        para <- list()
        k=0
        int=NULL
        n.int=NULL
        if (verbose==TRUE){
          cat("Checking for signaling between cell types",fill=TRUE)
        }
        for (i in z){
          if (sum(cluster==i)>1){
            tmp <- data[,cluster==i]
            tmp <- tmp[rowSums(tmp)>0,]
            if (sum(is.element(lig, rownames(tmp)))>0){
              lig.tmp <- rownames(tmp)[is.element(rownames(tmp),lig)]
              #lig.tmp <- lig.tmp[!is.element(lig.tmp,gene.list[[i]])]
            } else {lig.tmp=NULL}

            final.tmp <- SingleCellSignalR_LRDB[is.element(SingleCellSignalR_LRDB$ligand,lig.tmp),seq_len(2)]
            final.tmp <- data.frame(final.tmp,as.character(
              rep("paracrine",sum(is.element(SingleCellSignalR_LRDB$ligand,lig.tmp)))),
              stringsAsFactors=FALSE)
            m.lig <- rowSums(tmp[unique(final.tmp[,1]),])/sum(cluster==i)
            names(m.lig) <- unique(final.tmp[,1])

            if (sum(is.element(rec, rownames(tmp)))>0){
              rec.tmp <- rownames(tmp)[is.element(rownames(tmp[
                apply(tmp,1,function(x) sum(x>0))>tol*ncol(tmp),]),rec)]
            } else {rec.tmp=NULL}

            for (j in z[-i]){
              if (sum(cluster==j)>1){
                temp <- data[,cluster==j]
                temp <- temp[rowSums(temp)>0,]
                if (sum(is.element(rec, rownames(temp)))>0){
                  rec.temp <- rownames(temp)[is.element(rownames(temp),rec)]
                  #rec.temp <- rec.temp[!is.element(rec.temp,gene.list[[j]])]
                } else {rec.temp=NULL}
                rec.temp <- rec.temp[!is.element(rec.temp,rec.tmp)]
                m.rec <- rowSums(data.frame(temp[rec.temp,]))/sum(cluster==j)
                names(m.rec) <- rec.temp

                final <- final.tmp[is.element(final.tmp$receptor,rec.temp),]
                final <- cbind(final,SingleCellSignalR_LRscore(m.lig[final$ligand],m.rec[final$receptor],
                                                               med))
                exclus <- final$ligand %in% gene.list[[i]] & final$receptor %in%
                  gene.list[[j]]
                if (sum(exclus)!=0){
                  f.exclu <- final[exclus,]
                  final <- final[!(final$ligand %in% gene.list[[i]] & final$receptor
                                   %in% gene.list[[j]]),]
                  f.exclu[,3] <- "specific"
                  final <- rbind(f.exclu,final)
                }

                colnames(final) <- c(c.names[i],c.names[j],"interaction type",
                                     "SingleCellSignalR_LRscore")
                final <- final[final[,4]>s.score,]
                final <- final[order(final[,4],decreasing=TRUE),]

                if (species=="mus musculus"){
                  final[,1] <- Hs2mm[as.character(final[,1])]
                  final[,2] <- Hs2mm[as.character(final[,2])]
                }

                if (nrow(final)>0){
                  k=k+1
                  para[[k]] <- final
                  if (verbose==TRUE){
                    cat(paste(nrow(final),"interactions from",c.names[i],"to",
                              c.names[j]),fill=TRUE)
                  }
                  int <- c(int,paste(i,"-",j,sep=""))
                  n.int <- c(n.int,paste(c.names[i],"-",c.names[j],sep=""))
                  gr <- graph_from_data_frame(final,directed=FALSE)
                  if (write==TRUE){
                    fwrite(data.frame(final),paste(
                      "./cell-signaling/LR_interactions_",c.names[i],"-",c.names[j],
                      "-",int.type,".txt",sep=""),sep="\t")
                  }
                } else {
                  if (verbose==TRUE){
                    cat(paste(nrow(final),"No significant interaction found from",c.names[i],"to",
                              c.names[j]),fill=TRUE)
                  }
                }
              }
            }
          }
        }
        if (k!=0){
          names(para) <- n.int
        }
      }

      ## Returns ---------------------
      if (int.type=="autocrine"){
        return(auto)
      }
      if (int.type=="paracrine"){
        return(para)
      }
    }
    SingleCellSignalR_LRscore <- function(l,r,s){
      L=l^(1/2)
      R=r^(1/2)
      S=s
      sc=L*R/(S+L*R)
      return(sc)
    }
    SingleCellSignalR_LRDB <- LR_data
    SingleCellSignalR_LRDB$source <- c(rep('other',nrow(SingleCellSignalR_LRDB)))
    SingleCellSignalR_LRDB$PMIDs <-c(rep('other',nrow(SingleCellSignalR_LRDB)))
    colnames(SingleCellSignalR_LRDB) <- colnames(SingleCellSignalR::LRdb)
    signal <- SingleCellSignalR_cell_signaling(gene_expr_dataframe, rownames(gene_expr_dataframe), cluster = cluster1$cluster,c.names = NULL,species ="homo sapiens",tol=0, write = FALSE) # 信号分析执行细胞间网络分析
    inter.net <- inter_network(gene_expr_dataframe, rownames(gene_expr_dataframe), cluster = cluster1$cluster,signal,c.names = NULL, species ="homo sapiens",write = FALSE) # 细胞间网络分析
    rm(gene_expr_dataframe)
    SingleCellSignalR <- inter.net[["full-network"]]
    SingleCellSignalR_ligand_type <- as.data.frame(matrix(unlist(str_split(SingleCellSignalR$ligand,"\\.")),nrow = length(str_split(SingleCellSignalR$ligand,"\\.")),byrow = T))
    colnames(SingleCellSignalR_ligand_type) <- c("source","ligand")
    SingleCellSignalR_receptor_type <- as.data.frame(matrix(unlist(str_split(SingleCellSignalR$receptor,"\\.")),nrow = length(str_split(SingleCellSignalR$receptor,"\\.")),byrow = T))
    colnames(SingleCellSignalR_receptor_type) <- c("target","receptor")
    SingleCellSignalR1 <- data.frame(SingleCellSignalR_ligand_type,SingleCellSignalR_receptor_type,SingleCellSignalR[,"SingleCellSignalR_LRscore"])
    cell_type_cluster$cluster <- paste0("cluster"," ",cell_type_cluster$cluster)
    cell_type_cluster1 <- cell_type_cluster
    colnames(cell_type_cluster1)[2] <- c("source")
    cell_type_cluster2 <- cell_type_cluster
    colnames(cell_type_cluster2)[2] <- c("target")
    SingleCellSignalR2 <- merge(SingleCellSignalR1,cell_type_cluster1,by='source')
    SingleCellSignalR2 <- merge(SingleCellSignalR2,cell_type_cluster2,by='target')
    SingleCellSignalR3 <- SingleCellSignalR2[,c("cell_type_final.x","ligand","receptor","cell_type_final.y","SingleCellSignalR[, \"SingleCellSignalR_LRscore\"]")]
    colnames(SingleCellSignalR3) <- c("source","ligand","receptor","target","LRscore")
    SingleCellSignalR3$interaction_name <- paste0(SingleCellSignalR3$ligand,"_",SingleCellSignalR3$receptor)
    SingleCellSignalR3 <- merge(SingleCellSignalR3,db,by="interaction_name")
    SingleCellSignalR_net <- SingleCellSignalR3[,c("source","ligand.x","receptor.x","target","interaction_name")]
    colnames(SingleCellSignalR_net) <- c("source","ligand","receptor","target","interaction_name")
    print('SingleCellSignalR is over')
    #####Return
    ##Give a label to the result obtained by each inference method1
    cellchat_net$CellChat <- 1
    cellTalker_net$cellTalker <- 1
    Connectome_net$Connectome <- 1
    iTALK_net$iTALK <- 1
    ICELLNET_net$ICELLNET <- 1
    scSeqComm_net$scSeqComm <- 1
    SingleCellSignalR_net$SingleCellSignalR <- 1
    ##Results obtained by combining methods
    cellchat_cellTalker_net <- full_join(cellchat_net,cellTalker_net,by=c("source","ligand","receptor","target","interaction_name"),relationship = "many-to-many")
    cellchat_cellTalker_Connectome_net <- full_join(cellchat_cellTalker_net,Connectome_net,by=c("source","ligand","receptor","target","interaction_name"),relationship = "many-to-many")
    cellchat_cellTalker_Connectome_iTALK_net <- full_join(cellchat_cellTalker_Connectome_net,iTALK_net,by=c("source","ligand","receptor","target","interaction_name"),relationship = "many-to-many")
    cellchat_cellTalker_Connectome_iTALK_ICELLNET_net <- full_join(cellchat_cellTalker_Connectome_iTALK_net,ICELLNET_net,by=c("source","ligand","receptor","target","interaction_name"),relationship = "many-to-many")
    cellchat_cellTalker_Connectome_iTALK_ICELLNET_scSeqComm_net <- full_join(cellchat_cellTalker_Connectome_iTALK_ICELLNET_net,scSeqComm_net,by=c("source","ligand","receptor","target","interaction_name"),relationship = "many-to-many")
    cellchat_cellTalker_Connectome_iTALK_ICELLNET_scSeqComm_SingleCellSignalR_net <- full_join(cellchat_cellTalker_Connectome_iTALK_ICELLNET_scSeqComm_net,SingleCellSignalR_net,by=c("source","ligand","receptor","target","interaction_name"),relationship = "many-to-many")
    cellchat_cellTalker_Connectome_iTALK_ICELLNET_scSeqComm_SingleCellSignalR_net[is.na(cellchat_cellTalker_Connectome_iTALK_ICELLNET_scSeqComm_SingleCellSignalR_net)] <- 0
    CCC_net <- cellchat_cellTalker_Connectome_iTALK_ICELLNET_scSeqComm_SingleCellSignalR_net
    CCC_net$Freq <- rowSums(CCC_net[,6:ncol(CCC_net)])
    ##return
    return(CCC_net)
  }
  if(is.null(ligand_receptor_DB)==T){
    #####Run cellchat
    print('CellChat is running')
    cellchat <- createCellChat(object = sc_data, meta = sc_data@meta.data, group.by = "cell_type_final" )
    cellchat@DB <- CellChatDB.human
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat_net <- subsetCommunication(cellchat)
    cellchat_net <- cellchat_net[,c("source","ligand","receptor","target","interaction_name")]
    print('CellChat is over')
    #####Run scSeqComm
    print('scSeqComm is running')
    gene_expr_matrix <- as.matrix(sc_count)
    cell_metadata <- data.frame(colnames(as.matrix(sc_count)),sc_metadata$cell_type_final)
    colnames(cell_metadata) <- c("Cell_ID","Cluster_ID")
    scSeqComm_LR_DB <- scSeqComm::available_LR_pairs()[1:16]
    scSeqComm_LR_DB1 <- NULL
    for(i in 1:16){
      x=get(scSeqComm_LR_DB[i])[,c('ligand','receptor')]
      scSeqComm_LR_DB1 <- rbind(scSeqComm_LR_DB1,x)
      print(i)
    }
    scSeqComm_LR_DB1 <- scSeqComm_LR_DB1[!duplicated(scSeqComm_LR_DB1),]
    scSeqComm <- scSeqComm_analyze(gene_expr_matrix,
                                        cell_metadata = cell_metadata,
                                        inter_signaling = TRUE,
                                        LR_pairs_DB =scSeqComm_LR_DB1,
                                        intra_signaling = FALSE,
                                        N_cores = 1)
    rm(gene_expr_matrix)
    scSeqComm_result <- scSeqComm[['comm_results']]
    scSeqComm_net <- scSeqComm_result[,c("cluster_L","ligand","receptor","cluster_R")]
    scSeqComm_net$interaction_name <- paste0(scSeqComm_net$ligand,"_",scSeqComm_net$receptor)
    colnames(scSeqComm_net) <- c("source","ligand","receptor","target","interaction_name")
    print('scSeqComm is over')
    #####Run cellTalker
    print('cellTalker is running')
    cellTalker <- celltalk(input_object=sc_data,
                                metadata_grouping="cell_type_final",
                                ligand_receptor_pairs=ramilowski_pairs,
                                number_cells_required=100,
                                min_expression=1000,
                                max_expression=20000,
                                scramble_times=10)
    cellTalker_net <- cellTalker[,c("cell_type1","cell_type2","interaction")]
    cellTalker_net_LR <- matrix(unlist(strsplit(cellTalker_net$interaction,"_")),nrow = length(cellTalker_net$interaction),byrow = T)
    colnames(cellTalker_net_LR) <- c("ligand",'receptor')
    cellTalker_net_LR <- as.data.frame(cellTalker_net_LR)
    cellTalker_net <- cbind(cellTalker_net,cellTalker_net_LR)
    cellTalker_net <- cellTalker_net[,c("cell_type1","ligand",'receptor',"cell_type2","interaction")]
    colnames(cellTalker_net) <- c("source","ligand","receptor","target","interaction_name")
    print('cellTalker is over')
    #####Run Connectome
    print('Connectome is running')
    Connectome <- CreateConnectome(sc_data,
                                        LR.database = "fantom5",
                                        species = 'human')
    Connectome <- Connectome[-which(Connectome$recept.expression==0),]
    Connectome <- Connectome[-which(Connectome$ligand.expression==0),]
    Connectome_net <- Connectome[,c("source","ligand","receptor","target")]
    Connectome_net$interaction_name <- paste0(Connectome_net$ligand,"_",Connectome_net$receptor)
    print('Connectome is over')
    #####Run iTALK
    print('iTALK is running')
    gene_expr <- as.data.frame(t(as.matrix(sc_data@assays$RNA@counts)))
    gene_expr$cell_type <- sc_data$cell_type_final
    gene_expr$compare_group <- sc_data$cell_type_final
    highly_gene_exprrs_genes<-rawParse(gene_expr,top_genes=50,stats='mean')
    comm_list<-c('growth factor','other','cytokine','checkpoint')
    iTALK<-NULL
    for(i in comm_list){
      res_cat<-FindLR(highly_gene_exprrs_genes,datatype='mean count',comm_type=i)
      res_cat<-res_cat[order(res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs,decreasing=T),]
      iTALK<-rbind(iTALK,res_cat)
      print(i)
    }
    iTALK<-iTALK[order(iTALK$cell_from_mean_exprs*iTALK$cell_to_mean_exprs,decreasing=T),]
    rm(gene_expr)
    iTALK <- iTALK[iTALK[["cell_from_mean_exprs"]] != 0, ]
    iTALK <- iTALK[iTALK[["cell_to_mean_exprs"]] != 0, ]
    iTALK_net <- iTALK[,c("cell_from","ligand","receptor","cell_to")]
    colnames(iTALK_net) <- c("source","ligand","receptor","target")
    iTALK_net$interaction_name <- paste0(iTALK_net$ligand,"_",iTALK_net$receptor)
    print('iTALK is over')
    #####Run icellnet
    print('icellnet is running')
    ##load require LR_data
    data('ICELLNETdb',package = 'MLRS')
    icellnet_average_clean= sc.data.cleaning(object = sc_data, db = ICELLNETdb, filter.perc = 0, save_file = F, path="", force.file = F)
    icellnet_data=as.data.frame(gene.scaling(as.data.frame(icellnet_average_clean), n=1, db=ICELLNETdb))
    icellnet_PC.target=data.frame("Class"=unique(sc_metadata$cell_type_final), "ID"= unique(sc_metadata$cell_type_final), "Cell_type"=unique(sc_metadata$cell_type_final))
    rownames(icellnet_PC.target)=icellnet_PC.target$Cell_type
    icellnet <- NULL
    for(i in icellnet_PC.target$Cell_type){
      score.computation= icellnet.score(direction="out", PC.data=icellnet_data,
                                        CC.data= as.data.frame(icellnet_data[,i], row.names = rownames(icellnet_data)),
                                        PC.target = icellnet_PC.target, PC=icellnet_PC.target$Cell_type, CC.type = "RNAseq",
                                        PC.type = "RNAseq",db=ICELLNETdb)
      LR_score=data.frame(score.computation[[2]])
      colnames(LR_score) <- icellnet_PC.target$Cell_type
      LR_score$cell_type <- i
      icellnet <- rbind(icellnet,LR_score)
      score.computation <- c()
      print(i)
    }
    icellnet1 <- icellnet
    icellnet1$ligand <- ICELLNETdb$`Ligand 1`
    icellnet1$receptor <- ICELLNETdb$`Receptor 1`
    colnames(icellnet1) <- c(icellnet_PC.target$Cell_type,'source_cell_type',"ligand",'receptor')
    icellnet_net <- NULL
    for(i in icellnet_PC.target$Cell_type){
      x <- data.frame(icellnet1[,c("source_cell_type","ligand","receptor")],icellnet1[,i])
      x$target_cell_type <- i
      icellnet_net <- rbind(icellnet_net,x)
      x <- c()
      print(i)
    }
    colnames(icellnet_net) <- c("source_cell_type","ligand","receptor",'LR_score','target_cell_type')
    icellnet_net <- icellnet_net[,c("source_cell_type","ligand","receptor",'target_cell_type','LR_score')]
    icellnet_net <- na.omit(icellnet_net)
    icellnet_net <- icellnet_net[-which(icellnet_net$LR_score==0),]
    rownames(icellnet_net) <- 1:nrow(icellnet_net)
    icellnet_net$interaction_name <- paste0(icellnet_net$ligand,"_",icellnet_net$receptor)
    ICELLNET_net <- icellnet_net[,c("source_cell_type","ligand","receptor","target_cell_type","interaction_name")]
    colnames(ICELLNET_net) <- c("source","ligand","receptor","target","interaction_name")
    print('icellnet is over')
    #####Run SingleCellSignalR
    print('SingleCellSignalR is running')
    gene_expr_dataframe <- as.data.frame(sc_count)
    cluster <- data.frame(sc_data$cell_type_final)
    colnames(cluster) <- c('cell_type_final')
    cell_type_cluster <- data.frame(unique(sc_data$cell_type_final))
    colnames(cell_type_cluster) <- c('cell_type_final')
    cell_type_cluster$cluster <- 1:nrow(cell_type_cluster)
    cluster1 <- inner_join(cluster,cell_type_cluster,by = join_by(cell_type_final))
    cluster1$cluster <- as.numeric(cluster1$cluster)
    signal <- cell_signaling(gene_expr_dataframe, rownames(gene_expr_dataframe), cluster = cluster1$cluster,c.names = NULL,species ="homo sapiens",tol=0, write = FALSE) # 信号分析执行细胞间网络分析
    inter.net <- inter_network(gene_expr_dataframe, rownames(gene_expr_dataframe), cluster = cluster1$cluster,signal,c.names = NULL, species ="homo sapiens",write = FALSE) # 细胞间网络分析
    rm(gene_expr_dataframe)
    SingleCellSignalR <- inter.net[["full-network"]]
    SingleCellSignalR_ligand_type <- as.data.frame(matrix(unlist(str_split(SingleCellSignalR$ligand,"\\.")),nrow = length(str_split(SingleCellSignalR$ligand,"\\.")),byrow = T))
    colnames(SingleCellSignalR_ligand_type) <- c("source","ligand")
    SingleCellSignalR_receptor_type <- as.data.frame(matrix(unlist(str_split(SingleCellSignalR$receptor,"\\.")),nrow = length(str_split(SingleCellSignalR$receptor,"\\.")),byrow = T))
    colnames(SingleCellSignalR_receptor_type) <- c("target","receptor")
    SingleCellSignalR1 <- data.frame(SingleCellSignalR_ligand_type,SingleCellSignalR_receptor_type,SingleCellSignalR[,"LRscore"])
    cell_type_cluster$cluster <- paste0("cluster"," ",cell_type_cluster$cluster)
    cell_type_cluster1 <- cell_type_cluster
    colnames(cell_type_cluster1)[2] <- c("source")
    cell_type_cluster2 <- cell_type_cluster
    colnames(cell_type_cluster2)[2] <- c("target")
    SingleCellSignalR2 <- merge(SingleCellSignalR1,cell_type_cluster1,by='source')
    SingleCellSignalR2 <- merge(SingleCellSignalR2,cell_type_cluster2,by='target')
    SingleCellSignalR3 <- SingleCellSignalR2[,c("cell_type_final.x","ligand","receptor","cell_type_final.y","SingleCellSignalR[, \"LRscore\"]")]
    colnames(SingleCellSignalR3) <- c("source","ligand","receptor","target","LRscore")
    SingleCellSignalR3$interaction_name <- paste0(SingleCellSignalR3$ligand,"_",SingleCellSignalR3$receptor)
    SingleCellSignalR_net <- SingleCellSignalR3[,c("source","ligand","receptor","target","interaction_name")]
    print('SingleCellSignalR is over')
    #####Return
    ##Give a label to the result obtained by each inference method1
    cellchat_net$CellChat <- 1
    cellTalker_net$cellTalker <- 1
    Connectome_net$Connectome <- 1
    iTALK_net$iTALK <- 1
    ICELLNET_net$ICELLNET <- 1
    scSeqComm_net$scSeqComm <- 1
    SingleCellSignalR_net$SingleCellSignalR <- 1
    ##Results obtained by combining methods
    cellchat_cellTalker_net <- full_join(cellchat_net,cellTalker_net,by=c("source","ligand","receptor","target","interaction_name"),relationship = "many-to-many")
    cellchat_cellTalker_Connectome_net <- full_join(cellchat_cellTalker_net,Connectome_net,by=c("source","ligand","receptor","target","interaction_name"),relationship = "many-to-many")
    cellchat_cellTalker_Connectome_iTALK_net <- full_join(cellchat_cellTalker_Connectome_net,iTALK_net,by=c("source","ligand","receptor","target","interaction_name"),relationship = "many-to-many")
    cellchat_cellTalker_Connectome_iTALK_ICELLNET_net <- full_join(cellchat_cellTalker_Connectome_iTALK_net,ICELLNET_net,by=c("source","ligand","receptor","target","interaction_name"),relationship = "many-to-many")
    cellchat_cellTalker_Connectome_iTALK_ICELLNET_scSeqComm_net <- full_join(cellchat_cellTalker_Connectome_iTALK_ICELLNET_net,scSeqComm_net,by=c("source","ligand","receptor","target","interaction_name"),relationship = "many-to-many")
    cellchat_cellTalker_Connectome_iTALK_ICELLNET_scSeqComm_SingleCellSignalR_net <- full_join(cellchat_cellTalker_Connectome_iTALK_ICELLNET_scSeqComm_net,SingleCellSignalR_net,by=c("source","ligand","receptor","target","interaction_name"),relationship = "many-to-many")
    cellchat_cellTalker_Connectome_iTALK_ICELLNET_scSeqComm_SingleCellSignalR_net[is.na(cellchat_cellTalker_Connectome_iTALK_ICELLNET_scSeqComm_SingleCellSignalR_net)] <- 0
    CCC_net <- cellchat_cellTalker_Connectome_iTALK_ICELLNET_scSeqComm_SingleCellSignalR_net
    CCC_net$Freq <- rowSums(CCC_net[,6:ncol(CCC_net)])
    ##return
    return(CCC_net)
  }
}
