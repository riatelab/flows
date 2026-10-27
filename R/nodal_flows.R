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
#' res[1:5, 1:5]
nodal_flows <- function(mat) {
  diag(mat) <- 0
  res <- select_flows(mat = mat, method = "nfirst", global = FALSE, k = 1)
  res2 <- select_flows(mat = mat, method = "dominant", w = colSums(mat), k = 1)
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
#' @param leg_pos Position of the legend
#' @param add A boolean, if TRUE, add the layer to an existing plot.
#' @export
#' @importFrom mapsf mf_get_links mf_map mf_get_breaks
#' @importFrom maplegend leg_comp leg_draw
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
#' map_nodal_flows(
#'   mat = mat, x = UA,
#'   col_node = c("red", "orange", "yellow"),
#'   col_flow = "grey30",
#'   breaks = c(4, 100, 1000, 2500, 8655),
#'   lwd = c(1, 4, 8, 16), add = TRUE
#' )
#' mf_title("Dominant flows")
map_nodal_flows <- function(mat,
                            x,
                            inches = .15,
                            col_node = c("red", "orange", "yellow"),
                            breaks = "equal",
                            nbreaks = 4,
                            lwd = c(1, 5, 10, 20),
                            col_flow = "grey20",
                            leg_node = c(
                              "Dominant",
                              "Intermediate",
                              "Dominated",
                              "Size proportional\nto sum of inflows"
                            ),
                            leg_flow = "Flow intensity",
                            leg_pos = "topleft",
                            add = FALSE) {
  op <- par(lend = 1)
  on.exit(par(op), add = TRUE)
  diag(mat) <- 0
  w <- colSums(mat)
  dt <- data.frame(id = names(w), w)
  names(x)[1] <- "id"
  x <- merge(x, dt, by.x = names(x)[1], by.y = names(dt)[1], all.x = TRUE)
  res <- select_flows(mat = mat, method = "nfirst", global = FALSE, k = 1)
  res2 <- select_flows(mat = mat, method = "dominant", w = w, k = 1)
  mm <- mat * res * res2
  fdom <- cbind(i = row.names(mm), stack(as.data.frame(mm)))
  names(fdom) <- c("i", "fij", "j")
  fdom <- fdom[fdom$fij > 0, ]
  links <- mf_get_links(
    x = x,
    df = fdom,
    x_id = names(x)[1],
    df_id = c("i", "j")
  )
  x$node <- NA
  x[x$id %in% fdom$j & !x$id %in% fdom$i, "node"] <- "Dominant"
  x[x$id %in% fdom$j & x$id %in% fdom$i, "node"] <- "Intermediate"
  x[!x$id %in% fdom$j & x$id %in% fdom$i, "node"] <- "Dominated"
  x <- x[!is.na(x$node), ]



  mf_map(links, "fij", "grad",
    add = add, col = col_flow,
    nbreaks = nbreaks, breaks = breaks,
    leg_pos = NA, leg_title = leg_flow,
    lwd = lwd, leg_val_rnd = 0
  )
  mf_map(x, c("w", "node"), "prop_typo",
    add = TRUE, inches = inches, lwd = 1,
    pal = col_node, leg_pos = NA,
    val_order = c("Dominant", "Intermediate", "Dominated"),
    leg_val_rnd = c(0, 0)
  )

  breaks <- mf_get_breaks(x = links[["fij"]], nbreaks = nbreaks,
                          breaks = breaks)
  l1 <- list(
    type = "grad_line",
    val =  breaks,
    col = col_flow,
    lwd = lwd,
    title = leg_flow,
    val_rnd = 0
  )
  l1 <- do.call(leg_comp, l1)
  l2 <- list(
    leg = l1,
    type = "symb",
    pch = rep(21, 4),
    cex = c(2.8, 2, 1, 1),
    title = NA,
    val = leg_node,
    border = c(rep(getOption("mapsf.background"), 3), NA),
    lwd = 1,
    pal = c(col_node, NA),
    box_cex = c(1, .75)
  )
  l2 <- do.call(leg_comp, l2)
  leg_draw(l2, pos = leg_pos,mar = getOption("mapsf.mar"), frame = FALSE)

  return(invisible(list(nodes = x, flows = links)))
}




#' @title Nodal flows graph
#' @name plot_nodal_flows
#' @description This function plots a dominant flows graph.
#' @param mat A square matrix of dominant flows (see \link{nodal_flows}).
#' @param leg_pos Position of the legend, one of "topleft", "top",
#' "topright", "left", "right", "bottomleft", "bottom", "bottomright".
#' @param leg_flow Title of the flows legend.
#' @param leg_node Title of the nodes legend.
#' @param leg_node_modalities Text of node modalities.
#' @param labels A boolean, if TRUE, labels of dominant and intermediary nodes are plotted.
#' @note As square matrices can easily be plotted with plot.igraph or
#' gplot functions from igraph and sna packages respectively, we do not propose
#' visualisation for other output.
#' @details This function uses the Davidson Harel algorithm from igraph.
#' @returns The graph is returned, a graph is displayed.
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
                             leg_pos = "topright",
                             leg_flow = "Flow Intensity",
                             leg_node = "Node type",
                             leg_node_modalities = c("Dominant", "Intermediary", "Dominated"),
                             labels = FALSE) {
  g <- igraph::graph.adjacency(adjmatrix = mat, mode = "directed", weighted = TRUE)
  g <- igraph::delete.vertices(g, names(degree(g)[degree(g) == 0]))
  vertexdf <- data.frame(id = V(g)$name, col = NA, size = NA, name = NA)
  # Dominant
  vertexdf[(degree(g, mode = "in") > 0) & (degree(g, mode = "out") < 1), "col"] <- "red"
  vertexdf[(degree(g, mode = "in") > 0) & (degree(g, mode = "out") < 1), "size"] <- 2.8
  vertexdf[(degree(g, mode = "in") > 0) & (degree(g, mode = "out") < 1), "name"] <-
    as.character(vertexdf[(degree(g, mode = "in") > 0) & (degree(g, mode = "out") < 1), "id"])
  # intermediaire
  vertexdf[(degree(g, mode = "in") > 0) & (degree(g, mode = "out") > 0), "col"] <- "orange"
  vertexdf[(degree(g, mode = "in") > 0) & (degree(g, mode = "out") > 0), "size"] <- 2
  vertexdf[(degree(g, mode = "in") > 0) & (degree(g, mode = "out") > 0), "name"] <-
    as.character(vertexdf[(degree(g, mode = "in") > 0) & (degree(g, mode = "out") > 0), "id"])
  # Dominé
  vertexdf[(degree(g, mode = "in") < 1) & (degree(g, mode = "out") > 0), "col"] <- "yellow"
  vertexdf[(degree(g, mode = "in") < 1) & (degree(g, mode = "out") > 0), "size"] <- 1

  V(g)$color <- vertexdf$col
  V(g)$size <- vertexdf$size
  V(g)$names <- as.character(vertexdf$name)

  layout <-  igraph::layout.davidson.harel(g)
  # Normalize the layout in which the vertices are positioned
  vertex.layout <- edge.layout <- igraph::norm_coords(
    layout,
    ymin = -1,
    ymax = 1,
    xmin = -1,
    xmax = 1
  )

  E(g)$width <- ((E(g)$weight /  max(E(g)$weight))) * 10
  # Legend values and thickness
  breaks.edge <- seq(max(E(g)$weight), min(E(g)$weight), length.out = 4)


  plot.new()
  # Setting new graphic parameters
  opar <- par(usr = c(-1.1, 1.4, -1.1, 1.1), xpd = TRUE)
  on.exit(par(opar))

  # Display edges
  igraph::plot.igraph(
    g,
    edge.curved = .1,
    layout = edge.layout,
    vertex.size = 0,
    edge.color = "black",
    add = TRUE, edge.arrow.size = 0,
    edge.width = E(g)$width,
    vertex.label = NA
  )

  points(
    x = vertex.layout[, 1],
    y = vertex.layout[, 2],
    pch = 21,
    cex =  V(g)$size ,
    bg =  V(g)$color ,
    col = "grey",
    lwd = 1
  )

  if(isTRUE(labels)){
    text(
      x = edge.layout[, 1],
      y = edge.layout[, 2] + yinch(.01),
      labels = V(g)$names,
      family = "sans",
      col = "black",
      cex = .8,
      adj = c(.5, 0)
    )
  }



  l1 <- list(
    type = "prop_line",
    val = rev(breaks.edge),
    col = "black",
    lwd = max(E(g)$width),
    title = leg_flow,
    val_rnd = 0
  )
  l1 <- do.call(leg_comp, l1)
  l2 <- list(
    leg = l1,
    type = "symb",
    pch = rep(21, 3),
    cex = c(2.8, 2, 1),
    title = leg_node,
    val = leg_node_modalities,
    border = rep("black", 3),
    lwd = 1,
    pal = c("red", "orange", "yellow"),
    box_cex = c(1, .75)
  )
  l2 <- do.call(leg_comp, l2)
  leg_draw(l2, pos = leg_pos, mar = par("mar"), frame = FALSE)

  return(invisible(g))
}
