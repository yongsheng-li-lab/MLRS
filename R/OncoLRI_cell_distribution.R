#' @title Visualization of oncoLRI distribution in cells based on the CellChat package.
#'
#' @param OncoLRI Preferred oncoLRI containing at least the columns source, ligand, receptor, target, and the like.
#' @param color.use Colors represent different cell groups
#' @param title.name the name of the title
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param idents.use a vector giving the index or the name of cell groups of interest.
#' @param remove.isolate whether remove the isolate nodes in the communication network
#' @param top the fraction of interactions to show
#' @param weight.scale whether scale the weight
#' @param vertex.weight The weight of vertex: either a scale value or a vector
#' @param vertex.weight.max the maximum weight of vertex; defualt = max(vertex.weight)
#' @param vertex.size.max the maximum vertex size for visualization
#' @param vertex.label.cex The label size of vertex
#' @param vertex.label.color The color of label for vertex
#' @param edge.weight.max the maximum weight of edge; defualt = max(net)
#' @param edge.width.max The maximum edge width for visualization
#' @param alpha.edge the transprency of edge
#' @param label.edge Whether or not shows the label of edges
#' @param edge.label.color The color for single arrow
#' @param edge.label.cex The size of label for arrows
#' @param edge.curved Specifies whether to draw curved edges, or not. This can be a logical or a numeric vector or scalar. First the vector is replicated to have the same length as the number of edges in the graph. Then it is interpreted for each edge separately. A numeric value specifies the curvature of the edge; zero curvature means straight edges, negative values means the edge bends clockwise, positive values the opposite. TRUE means curvature 0.5, FALSE means curvature zero
#' @param shape The shape of the vertex, currently “circle”, “square”, “csquare”, “rectangle”, “crectangle”, “vrectangle”, “pie” (see vertex.shape.pie), ‘sphere’, and “none” are supported, and only by the plot.igraph command. “none” does not draw the vertices at all, although vertex label are plotted (if given). See shapes for details about vertex shapes and vertex.shape.pie for using pie charts as vertices.
#' @param layout The layout specification. It must be a call to a layout specification function.
#' @param margin The amount of empty space below, over, at the left and right of the plot, it is a numeric vector of length four. Usually values between 0 and 0.5 are meaningful, but negative values are also possible, that will make the plot zoom in to a part of the graph. If it is shorter than four then it is recycled.
#' @param vertex.size Deprecated. Use 'vertex.weight'
#' @param arrow.width The width of arrows
#' @param arrow.size the size of arrow
#' @param text.x the x-coordinates to add the text
#' @param text.y the y-coordinates to add the text
#'
#' @importFrom CellChat netVisual_circle
#'
#' @return
#' @export
#'
#' @examples
OncoLRI_cell_distribution <- function(OncoLRI,
                                      color.use = NULL,
                                      title.name = NULL,
                                      sources.use = NULL,
                                      targets.use = NULL,
                                      idents.use = NULL,
                                      remove.isolate = FALSE,
                                      top = 1,
                                      weight.scale = FALSE,
                                      vertex.weight = 20,
                                      vertex.weight.max = NULL,
                                      vertex.size.max = NULL,
                                      vertex.label.cex = 1,
                                      vertex.label.color = "black",
                                      edge.weight.max = NULL,
                                      edge.width.max = 8,
                                      alpha.edge = 0.6,
                                      label.edge = FALSE,
                                      edge.label.color = "black",
                                      edge.label.cex = 0.8,
                                      edge.curved = 0.2,
                                      shape = "circle",
                                      layout = in_circle(),
                                      margin = 0.2,
                                      vertex.size = NULL,
                                      arrow.width = 1,
                                      arrow.size = 0.2,
                                      text.x = 0,
                                      text.y = 1.5){
  if(isTRUE(unique(OncoLRI$source)=="Tumor")==T){
    OncoLRI_type <- rownames(table(OncoLRI[which(OncoLRI$source=="Tumor"),"target"]))
    OncoLRI_count <- matrix(nrow = length(unique(OncoLRI$target)),ncol = length(unique(OncoLRI$target)))
    dimnames(OncoLRI_count) <- list(OncoLRI_type,OncoLRI_type)
    OncoLRI_count['Tumor',] <- table(OncoLRI$target)
    OncoLRI_count[is.na(OncoLRI_count)] <- 0
    CellChat::netVisual_circle(OncoLRI_count,
                               color.use = color.use,
                               title.name = title.name,
                               sources.use = sources.use,
                               targets.use = targets.use,
                               idents.use = idents.use,
                               remove.isolate = remove.isolate,
                               top = top,
                               weight.scale = weight.scale,
                               vertex.weight = vertex.weight,
                               vertex.weight.max = vertex.weight.max,
                               vertex.size.max = vertex.size.max,
                               vertex.label.cex = vertex.label.cex,
                               vertex.label.color = vertex.label.color,
                               edge.weight.max = edge.weight.max,
                               edge.width.max = edge.width.max,
                               alpha.edge = alpha.edge,
                               label.edge = label.edge,
                               edge.label.color = edge.label.color,
                               edge.label.cex = edge.label.cex,
                               edge.curved = edge.curved,
                               shape = shape,
                               layout = layout,
                               margin = margin,
                               vertex.size = vertex.size,
                               arrow.width = arrow.width,
                               arrow.size = arrow.size,
                               text.x = text.x,
                               text.y = text.y)
  }
  if(isTRUE(unique(OncoLRI$target)=="Tumor")==T){
    OncoLRI_type <- rownames(table(OncoLRI[which(OncoLRI$target=="Tumor"),"source"]))
    OncoLRI_count <- matrix(nrow = length(unique(OncoLRI$source)),ncol = length(unique(OncoLRI$source)))
    dimnames(OncoLRI_count) <- list(OncoLRI_type,OncoLRI_type)
    OncoLRI_count[,'Tumor'] <- table(OncoLRI$source)
    OncoLRI_count[is.na(OncoLRI_count)] <- 0
    CellChat::netVisual_circle(OncoLRI_count,
                               color.use = color.use,
                               title.name = title.name,
                               sources.use = sources.use,
                               targets.use = targets.use,
                               idents.use = idents.use,
                               remove.isolate = remove.isolate,
                               top = top,
                               weight.scale = weight.scale,
                               vertex.weight = vertex.weight,
                               vertex.weight.max = vertex.weight.max,
                               vertex.size.max = vertex.size.max,
                               vertex.label.cex = vertex.label.cex,
                               vertex.label.color = vertex.label.color,
                               edge.weight.max = edge.weight.max,
                               edge.width.max = edge.width.max,
                               alpha.edge = alpha.edge,
                               label.edge = label.edge,
                               edge.label.color = edge.label.color,
                               edge.label.cex = edge.label.cex,
                               edge.curved = edge.curved,
                               shape = shape,
                               layout = layout,
                               margin = margin,
                               vertex.size = vertex.size,
                               arrow.width = arrow.width,
                               arrow.size = arrow.size,
                               text.x = text.x,
                               text.y = text.y)
  }
}
