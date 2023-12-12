#' @title Nodal flows selection
#' @name nodal_flows
#' @description Perform a Nystuen & Dacey's dominants flows analysis.
#' @param mat A square matrix of flows.
#' @export
#' @return
#' The matrix of the selected flows is returned.
#' @examples
#' nav <- read.csv(system.file("csv/nav.csv", package = "flows"))
#' mat <- prepare_mat(x = nav, i = "i", j = "j", fij = "fij")
#' res <- nodal_flows(mat)
#' res[1:5,1:5]
nodal_flows <- function(mat){
  diag(mat) <- 0
  res <- select_flows(mat = mat, method = "nfirst", global = FALSE, k = 1)
  res2 <- select_flows(mat = mat, method = "dominant", w =  colSums(mat), k = 1)
  mat * res * res2
}






#' @title Nodal flows map
#' @name map_nodal_flows
#' @description Perform a Nystuen & Dacey's dominants, or nodal, flows analysis and plot
#' a dominant flows map.
#' @param mat A square matrix of flows.
#' @param x An sf object, the first column contains a unique identifier
#' matching mat column and row names.
#' @param inches Size of the largest circle.
#' @param col_node Node colors, a vector of 3 colors.
#' @param breaks How to classify flows, either a numeric vector with the actual
#' breaks, or a classification method name (see mf_get_breaks())
#' @param nbreaks Number of classes.
#' @param lwd Flows widths
#' @param col_flow Flows color
#' @param leg_node Labels for the nodes legend
#' @param leg_flow Label for the flows legend
#' @param leg_pos_flow Position of the flows legend
#' @param leg_pos_node Position of the node legend
#' @param add A boolean, if TRUE, add the layer to an existing plot.
#' @export
#' @importFrom mapsf mf_get_links mf_map
#' @importFrom utils stack
#' @return A list of sf objects is returned. The first element contains the
#' nodes with their weight and classification (dominant, intermediary, dominated).
#' The second element contains the flows (i, j, fij)
#' @examples
#' library(sf)
#' library(mapsf)
#' nav <- read.csv(system.file("csv/nav.csv", package = "flows"))
#' mat <- prepare_mat(x = nav, i = "i", j = "j", fij = "fij")
#' UA <- st_read(system.file("gpkg/GE.gpkg", package = "flows"), layer = "urban_area")
#' GE <- st_read(system.file("gpkg/GE.gpkg", package = "flows"), layer = "region")
#' mf_map(GE)
#' map_nodal_flows(mat = mat, x = UA,
#'                  col_node = c('red', 'orange', 'yellow'),
#'                  col_flow = "grey30",
#'                  breaks = c(4,100,1000,2500,8655),
#'                  lwd = c(1,4,8,16), add = TRUE)
#' mf_title("Dominant flows")
map_nodal_flows <- function(mat,
                            x,
                            inches = .15,
                            col_node = c("red", "orange", "yellow"),
                            breaks = "equal",
                            nbreaks = 4,
                            lwd = c(1,5,10,20),
                            col_flow = "grey20",
                            leg_node = c("Dominant",
                                         "Intermediate",
                                         "Dominated",
                                         "Size proportional\nto sum of inflows"),
                            leg_flow = "Flow intensity",
                            leg_pos_flow = "topleft",
                            leg_pos_node = "topright",
                            add = FALSE){

  op <- par(lend = 1)
  on.exit(par(op), add = TRUE)
  diag(mat) <- 0
  w <- colSums(mat)
  dt <- data.frame(id = names(w), w)
  x <- merge(x, dt, by.x = names(x)[1], by.y = names(dt)[1], all.x = T)
  res <- select_flows(mat = mat, method = "nfirst", global = FALSE, k = 1)
  res2 <- select_flows(mat = mat, method = "dominant", w = w, k = 1)
  mm <- mat * res * res2
  fdom <- cbind(i = row.names(mm), stack(as.data.frame(mm)))
  names(fdom) <- c("i", "fij", "j")
  fdom <- fdom[fdom$fij > 0,]
  links <- mf_get_links(x = x, df = fdom, x_id = names(x)[1],
                        df_id = c("i", "j"))
  x$node <- NA
  x[x$id %in% fdom$j & !x$id %in% fdom$i, "node"] <- "Dominant"
  x[x$id %in% fdom$j & x$id %in% fdom$i, "node"] <- "Intermediate"
  x[!x$id %in% fdom$j & x$id %in% fdom$i, "node"] <- "Dominated"
  x <- x[x$node != "green",]

  mf_map(links, "fij", "grad", add = add, col = col_flow,
         nbreaks = nbreaks, breaks = breaks,
         leg_pos = leg_pos_flow, leg_title = leg_flow,
         lwd = lwd, leg_val_rnd = 0)
  mf_map(x, c("w", "node"), "prop_typo", add = T, inches = inches,
         pal = col_node, leg_pos = c(NA,NA),
         val_order = c("Dominant", "Intermediate", "Dominated"),
         leg_val_rnd = c(0,0)
  )
  # print(leg_pos_node)
  mapsf::mf_legend(type = "symb", pch = rep(21,4),
               cex = c(2.8,2,1,1), title = NA,
               val = leg_node,
               pos = leg_pos_node, border = c(1,1,1,NA),
               pal = c(col_node, NA), box_cex = c(1,.5))



  return(invisible(list(nodes = x, flows = links)))

}




#' @title Nodal flows graph
#' @name plot_nodal_flows
#' @description This function plots a dominant flows graph.
#' @param mat A square matrix of dominant flows (see \link{nodal_flows}).
#' @param leg_pos_flows Position of the flows legend, one of "topleft", "top",
#' "topright", "left", "right", "bottomleft", "bottom", "bottomright".
#' @param leg_flow Title of the flows legend.
#' @param leg_pos_node Position of the nodes legend, one of "topleft", "top",
#' "topright", "left", "right", "bottomleft", "bottom", "bottomright".
#' @param leg_node Text of the nodes legend.
#' @param labels A boolean, if TRUE, labels of dominant and intermediary nodes are plotted.
#' @note As square matrices can easily be plotted with \link[igraph]{plot.igraph} or
#' \link[sna]{gplot} functions from igraph and sna packages, we do not propose
#' visualisation for other outputs.
#' @details This function uses the Davidson Harel algorithm from igraph.
#' @import graphics
#' @importFrom igraph "V" "E" "V<-" "E<-" "degree"
#' @examples
#' nav <- read.csv(system.file("csv/nav.csv", package = "flows"))
#' mat <- prepare_mat(x = nav, i = "i", j = "j", fij = "fij")
#' res <- nodal_flows(mat)
#'
#' # Plot dominant flows graph
#' plot_nodal_flows(mat = res)
#' @export
plot_nodal_flows <- function(mat,
                             leg_pos_flows = "topright",
                             leg_flow = "Flows Intensity",
                             leg_pos_node = "bottomright",
                             leg_node = c("Dominant", "Intermediary",
                                          "Dominated",
                                          "Size proportional\nto sum of inflows"),
                             labels = FALSE){
  g <- igraph::graph.adjacency(adjmatrix = mat,mode = "directed", weighted = TRUE)
  g <- igraph::delete.vertices(g, names(degree(g)[degree(g)==0]))
  vertexdf <-  data.frame(id = V(g)$name, col = NA, size = NA, name = NA)
  # Dominant
  vertexdf[(degree(g, mode = "in") > 0) & (degree(g, mode = "out") < 1), "col"] <- "red"
  vertexdf[(degree(g, mode = "in") > 0) & (degree(g, mode = "out") < 1), "size"] <- 6
  vertexdf[(degree(g, mode = "in") > 0) & (degree(g, mode = "out") < 1), "name"] <-
    as.character(vertexdf[(degree(g, mode = "in") > 0) & (degree(g, mode = "out") < 1), "id"])
  # intermediaire
  vertexdf[(degree(g, mode = "in") > 0) & (degree(g, mode = "out") > 0), "col"] <- "orange"
  vertexdf[(degree(g, mode = "in") > 0) & (degree(g, mode = "out") > 0), "size"] <- 4
  vertexdf[(degree(g, mode = "in") > 0) & (degree(g, mode = "out") > 0), "name"]<-
    as.character(vertexdf[(degree(g, mode = "in") > 0) & (degree(g, mode = "out") > 0), "id"])
  # Dominé
  vertexdf[(degree(g, mode = "in") < 1) & (degree(g, mode = "out") > 0), "col"] <- "yellow"
  vertexdf[(degree(g, mode = "in") < 1) & (degree(g, mode = "out") > 0), "size"] <- 2

  V(g)$color <- vertexdf$col
  V(g)$size <- vertexdf$size
  V(g)$names <- as.character(vertexdf$name)
  E(g)$color <- "black"
  E(g)$width <- ((E(g)$weight) * 8 / (max(E(g)$weight)-min(E(g)$weight)))+1
  lg <- igraph::layout.davidson.harel(g)
  g <- igraph::set.graph.attribute(graph = g, name = "layout", value = lg)
  if(labels == TRUE){
    x <- igraph::plot.igraph(g, vertex.label = V(g)$names, vertex.label.cex = 1,
                             vertex.label.color = "black",
                             vertex.size = V(g)$size, edge.arrow.size = 0)
  }else{
    x <- igraph::plot.igraph(g, vertex.label = NA,
                             vertex.size = V(g)$size, edge.arrow.size = 0)
  }
  LegendPropLines(pos = leg_pos_flows, legTitle = leg_flow,
                  legTitleCex = 0.8, legValuesCex = 0.6,
                  varvect = c(min(E(g)$weight),max(E(g)$weight)),
                  sizevect = c(2, 10), col = "black",
                  frame = FALSE, round = 0)


  legend(x = leg_pos_node, legend = leg_node,
         cex = c(0.8), pt.cex = c(2.8,2,1,0), bty = "n",
         pt.bg = c("red", "orange", "yellow", NA),
         pch = c(21,21,21,21))

}
