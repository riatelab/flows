#' @title Package description
#' @name flows
#' @description Selections on flow matrices, statistics on selected flows, map
#' and graph visualisations.
#'
#' An introduction to the package conceptual background and usage is proposed in
#' a vignette (see \code{vignette(topic = "flows")}) and a paper (Beauguitte, Giraud & Guérois 2015).
#' @references L. Beauguitte, T. Giraud & M. Guérois, 2015. "Un outil pour la sélection et la visualisation de flux : le package flows",
#'  \emph{Netcom}, 29-3/4:399-408. \url{https://journals.openedition.org/netcom/2134}.
#' @docType package
NULL

#' @title Commuters datasets
#' @name commuters_datasets
#' @description
#' Data on commuters between Urban Areas of the French Grand Est region in 2011.
#' Fields: \cr
#' \itemize{
#' \item{i: Code of the urban area of residence}
#' \item{namei: Name of the urban area of residence}
#' \item{wi: Total number of active occupied persons in the urban area of residence}
#' \item{j: Code of the urban area of work}
#' \item{namej: Name of the urban area of work}
#' \item{wj: Total number of active occupied persons in the urban area of work}
#' \item{fij: Number of commuters between i and j}
#' }
#'
#' Geopackage of the Grand Est region in France and its urban areas (2010 delineation).
#'
#'
#' @references
#' Commuters dataset: \url{https://www.insee.fr/fr/statistiques/2022113} \cr
#' Spatial dataset: \url{https://www.data.gouv.fr/en/datasets/geofla-r}
#' @docType data
#' @examples
#' nav <- read.csv(system.file("csv/nav.csv", package = "flows"))
#' library(sf)
#' UA <- st_read(system.file("gpkg/GE.gpkg", package = "flows"), layer = "urban_area")
#' GE <- st_read(system.file("gpkg/GE.gpkg", package = "flows"), layer = "region")
NULL
