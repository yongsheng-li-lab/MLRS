######define two function
# function to shift the x-axis coordinates when points are too close
shift.lollipop.x <- function(mut.pos = NULL, total.length = NULL, shift.factor = 0.01){

  pos.dif <- 0
  for (i in 1:length(mut.pos)){
    pos.dif <- c(pos.dif, mut.pos[i+1] - mut.pos[i])
  }

  idx <- which(pos.dif < shift.factor*total.length)

  ## deal with odd and even sets of points
  if (median(idx) %% 1==0){
    mut.pos[idx[idx < median(idx)]] <- mut.pos[idx[idx < median(idx)]] - shift.factor*total.length
    mut.pos[idx[idx > median(idx)]] <- mut.pos[idx[idx > median(idx)]] + shift.factor*total.length
  } else {
    mut.pos[idx[idx == median(idx)-0.5]] <- mut.pos[idx[idx == median(idx)-0.5]] - 0.5*shift.factor*total.length
    mut.pos[idx[idx == median(idx)+0.5]] <- mut.pos[idx[idx == median(idx)+0.5]] + 0.5*shift.factor*total.length

    mut.pos[idx[idx < median(idx)-0.5]] <- mut.pos[idx[idx < median(idx)-0.5]] - shift.factor*total.length
    mut.pos[idx[idx > median(idx)+0.5]] <- mut.pos[idx[idx > median(idx)+0.5]] + shift.factor*total.length
  }

  mut.pos

}

# function to split the segment into 3 parts
shift.lollipop.y <- function(x, start.y = 0.7){
  mod.start <- x - start.y

  set1 <- start.y + mod.start/3
  set2 <- set1 + mod.start/3

  as.data.frame(cbind(set1,set2))
}
#####lollipop.Plot
#' @title Draws lollipop plot of interface residues on  Protein structure.
#'
#' @param gene_mutation Gene mutation data, including 4 columns, the first column AA indicates the gene amino acid site, the second column Mut is set to HGVSp_Short, the third column Type defaults to the mutation type, and the fourth column is the gene mutation count.
#' @param gene_residues_data Gene interface data, including 4 columns.
#' @param protein.df Protein structure length data.
#' @param protein_height map parameter
#' @param protein_fill map parameter
#' @param protein_color map parameter
#' @param protein_size map parameter
#' @param str.fill map parameter
#' @param str.col map parameter
#' @param residule.fill map parameter
#' @param residule.col map parameter
#' @param line.col map parameter
#' @param residule_height map parameter
#' @param linewidth map parameter
#' @param point_shape map parameter
#' @param point_cols map parameter
#' @param point_fill map parameter
#' @param point_stroke map parameter
#' @param point_size map parameter
#' @param is.lable Logical value, whether to add a label
#' @param str.size map parameter
#' @param residule.size map paramete
#' @param shift.factor map paramete
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggrepel geom_text_repel
#'
#' @return
#' @export
#'
#' @examples
lollipop.Plot <- function(gene_mutation,
                          gene_residues_data,
                          protein.df,
                          protein_height=0.1,
                          protein_fill=NA,
                          protein_color="#6b8e23",
                          protein_size=1,
                          str.fill="#cdbde2",
                          str.col="#cdbde2",
                          str.size=0,
                          residule.fill="#FF0000",
                          residule.col="#FF0000",
                          residule.size=0,
                          line.col="grey50",
                          residule_height=0.5,
                          linewidth= 0.5,
                          point_shape = 21,
                          point_cols,
                          point_fill,
                          point_stroke=0.2,
                          point_size = 3,
                          is.lable=F,
                          shift.factor=0.02
                         ){
  gene_mutation$Shift.AA <- shift.lollipop.x(gene_mutation$AA, range(gene_mutation$AA)[2],shift.factor = shift.factor)
  gene_mutation1 <- cbind(gene_mutation, shift.lollipop.y(gene_mutation$Freq, 0.7))
  if(is.lable==F){
    p <- ggplot() +
      geom_rect(data = subset(gene_residues_data, Type == "str"),
                mapping = aes(xmin =1, xmax = max(gene_mutation1$Shift.AA)+10),
                ymin = -residule_height, ymax = 0,
                fill = str.fill,
                size=str.size,
                colour = str.col)+
      scale_y_continuous(limits = c(-(residule_height+0.1),2))+
      scale_x_continuous(expand = c(0,0),limits = c(-floor((max(gene_mutation1$Shift.AA)+10)/6),max(gene_mutation1$Shift.AA)+10)
      )+
      geom_segment(data = gene_mutation1, linewidth=linewidth,colour=line.col,
                   mapping = aes(x = AA, xend = AA, y = 0, yend = 0.4)) +
      geom_segment(data = gene_mutation1,linewidth=linewidth,colour=line.col,
                   mapping = aes(x = AA, xend = Shift.AA, y = 0.4, yend = 0.8)) +
      geom_segment(data = gene_mutation1,linewidth=linewidth,colour=line.col,
                   mapping = aes(x = Shift.AA, xend = Shift.AA, y = 0.8, yend = 1)) +
      geom_rect(data = protein.df,
                mapping = aes(xmin =Start, xmax =End),
                ymin = -residule_height-0.1, ymax = protein_height,fill=protein_fill,colour = protein_color,size=protein_size)+
      geom_point(data = subset(gene_mutation, Type == names(point_cols)[2]),
                 mapping = aes(x = Shift.AA,size=Freq),
                 fill = point_fill[2],
                 colour = point_fill[2],
                 y = 1,
                 shape = point_shape,stroke=point_stroke
      )+
      geom_point(data = subset(gene_mutation, Type == names(point_cols)[1]),
                 mapping = aes(x = Shift.AA,size=Freq),
                 fill = point_fill[1],
                 colour = point_fill[1],
                 y = 1,
                 shape = point_shape,stroke=point_stroke
      ) +
      scale_size_continuous(range = c(point_size,point_size+max(gene_mutation$Freq)))+
      geom_rect(data = subset(gene_residues_data, Type == "residue"),
                mapping = aes(xmin = Start, xmax = End, ymin = -residule_height, ymax = 0, fill = Feature, group = Feature),
                fill = residule.fill[subset(gene_residues_data, Type == "residue")$Feature],
                size=residule.size,
                colour = residule.col)+
      labs(x=NULL,title = NULL,y=NULL)+
      theme_bw()+theme(panel.border = element_blank(),
                       panel.grid=element_blank(),
                       plot.margin = unit(c(0,0.5,0,0),'cm'),
                       legend.position = "none",
                       axis.title = element_blank(),
                       axis.text= element_blank(),
                       axis.ticks =element_blank(),
                       axis.ticks.length.y = unit(0.2,'cm'))
    print(p)
  }else{
    p <- ggplot() +
      geom_rect(data = subset(gene_residues_data, Type == "str"),
                mapping = aes(xmin =1, xmax = max(gene_mutation1$Shift.AA)+10),
                ymin = -residule_height, ymax = 0,
                fill = str.fill,
                size=str.size,
                colour = str.col)+
      scale_y_continuous(limits = c(-(residule_height*2),8))+
      scale_x_continuous(expand = c(0,0),limits = c(-floor((max(gene_mutation1$Shift.AA)+10)/6),max(gene_mutation1$Shift.AA)+10)
      )+
      geom_segment(data = gene_mutation1, linewidth=linewidth,colour=line.col,
                   mapping = aes(x = AA, xend = AA, y = 0, yend = 0.4)) +
      geom_segment(data = gene_mutation1,linewidth=linewidth,colour=line.col,
                   mapping = aes(x = AA, xend = Shift.AA, y = 0.4, yend = 0.8)) +
      geom_segment(data = gene_mutation1,linewidth=linewidth,colour=line.col,
                   mapping = aes(x = Shift.AA, xend = Shift.AA, y = 0.8, yend = 1)) +
      geom_rect(data = protein.df,
                mapping = aes(xmin =Start, xmax =End),
                ymin = -residule_height-0.1, ymax = protein_height,fill=protein_fill,colour = protein_color,size=protein_size)+
      geom_point(data = subset(gene_mutation, Type == names(point_cols)[2]),
                 mapping = aes(x = Shift.AA,size=Freq),
                 fill = point_fill[2],
                 colour = point_fill[2],
                 y = 1,
                 shape = point_shape,stroke=point_stroke
      )+
      geom_point(data = subset(gene_mutation, Type == names(point_cols)[1]),
                 mapping = aes(x = Shift.AA,size=Freq),
                 fill = point_fill[1],
                 colour = point_fill[1],
                 y = 1,
                 shape = point_shape,stroke=point_stroke
      ) +
      scale_size_continuous(range = c(point_size,point_size+max(gene_mutation$Freq)))+
      geom_rect(data = subset(gene_residues_data, Type == "residue"),
                mapping = aes(xmin = Start, xmax = End, ymin = -residule_height, ymax = 0, fill = Feature, group = Feature),
                fill = residule.fill[subset(gene_residues_data, Type == "residue")$Feature],
                size=residule.size,
                colour = residule.col)+
      labs(x=NULL,title = NULL,y=NULL)+
      geom_text_repel(data = gene_mutation1,
                      mapping = aes(x = Shift.AA, y = Freq, label = Mut),
                      bg.colour = "white",
                      seed = 12345,
                      nudge_y = 0.25,
                      angle = 90)+
      theme_bw()+theme(panel.border = element_blank(),
                       panel.grid=element_blank(),
                       plot.margin = unit(c(0,0.5,0,0),'cm'),
                       legend.position = "none",
                       axis.title = element_blank(),
                       axis.text= element_blank(),
                       axis.ticks =element_blank(),
                       axis.ticks.length.y = unit(0.2,'cm'))
    print(p)
  }
}
