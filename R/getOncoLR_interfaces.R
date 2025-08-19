#' @title Show which oncoLR occurs interfaces.
#'
#' @param OncoLR Preferred oncoLR.
#' @param mutation_data Mutation data.
#'
#' @return
#' @export
#'
#' @examples
getOncoLR_interfaces <- function(OncoLR,mutation_data){
  mutation_data1 <- mutation_data
  colnames(mutation_data1) <- c('ligand','gene_mutation_type','gene_interface_residues')
  ligand <- data.frame(OncoLR$ligand)
  colnames(ligand) <- c('ligand')
  ligand_mutation_interface <- merge(ligand,mutation_data1,by='ligand')
  ligand_mutation_interface$gene_type <- 'ligand'
  ligand_mutation_interface <- ligand_mutation_interface [,c('ligand','gene_type','gene_mutation_type','gene_interface_residues')]
  colnames(ligand_mutation_interface) <- c('gene','gene_type','gene_mutation_type','gene_interface_residues')
  mutation_data2 <- mutation_data
  colnames(mutation_data2) <- c('receptor','gene_mutation_type','gene_interface_residues')
  receptor <- data.frame(OncoLR$receptor)
  colnames(receptor) <- c('receptor')
  receptor_mutation_interface <- merge(receptor,mutation_data2,by='receptor')
  receptor_mutation_interface$gene_type <- 'receptor'
  receptor_mutation_interface <- receptor_mutation_interface [,c('receptor','gene_type','gene_mutation_type','gene_interface_residues')]
  colnames(receptor_mutation_interface) <- c('gene','gene_type','gene_mutation_type','gene_interface_residues')
  OncoLR_mutation_interface <- rbind(ligand_mutation_interface,receptor_mutation_interface)
  OncoLR_mutation_interface <- OncoLR_mutation_interface[!duplicated(OncoLR_mutation_interface),]
  return(OncoLR_mutation_interface)
}
