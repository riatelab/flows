% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nodal_flows.R
\name{map_nodal_flows}
\alias{map_nodal_flows}
\title{Nodal flows map}
\usage{
map_nodal_flows(
  mat,
  x,
  inches = 0.15,
  col_node = c("red", "orange", "yellow"),
  breaks = "equal",
  nbreaks = 4,
  lwd = c(1, 5, 10, 20),
  col_flow = "grey20",
  leg_node = c("Dominant", "Intermediate", "Dominated",
    "Size proportional\\nto sum of inflows"),
  leg_flow = "Flow intensity",
  leg_pos_flow = "topleft",
  leg_pos_node = "topright",
  add = FALSE
)
}
\arguments{
\item{mat}{A square matrix of flows.}

\item{x}{An sf object, the first column contains a unique identifier
matching mat column and row names.}

\item{inches}{Size of the largest circle.}

\item{col_node}{Node colors, a vector of 3 colors.}

\item{breaks}{How to classify flows, either a numeric vector with the actual
breaks, or a classification method name (see mf_get_breaks())}

\item{nbreaks}{Number of classes.}

\item{lwd}{Flows widths}

\item{col_flow}{Flows color}

\item{leg_node}{Labels for the nodes legend}

\item{leg_flow}{Label for the flows legend}

\item{leg_pos_flow}{Position of the flows legend}

\item{leg_pos_node}{Position of the node legend}

\item{add}{A boolean, if TRUE, add the layer to an existing plot.}
}
\value{
A list of sf objects is returned. The first element contains the
nodes with their weight and classification (dominant, intermediary, dominated).
The second element contains the flows (i, j, fij)
}
\description{
Perform a Nystuen & Dacey's dominants, or nodal, flows analysis and plot
a dominant flows map.
}
\examples{
library(sf)
library(mapsf)
nav <- read.csv(system.file("csv/nav.csv", package = "flows"))
mat <- prepare_mat(x = nav, i = "i", j = "j", fij = "fij")
UA <- st_read(system.file("gpkg/GE.gpkg", package = "flows"), layer = "urban_area")
GE <- st_read(system.file("gpkg/GE.gpkg", package = "flows"), layer = "region")
mf_map(GE)
map_nodal_flows(mat = mat, x = UA,
                 col_node = c('red', 'orange', 'yellow'),
                 col_flow = "grey30",
                 breaks = c(4,100,1000,2500,8655),
                 lwd = c(1,4,8,16), add = TRUE)
mf_title("Dominant flows")
}
