#' @title Flow matrix preparation
#' @name prepare_mat
#' @description From a long format matrix to a a wide format matrix.
#' @param x A data.frame of flows between origins and destinations: long format
#' matrix (origins, destinations, flows intensity).
#' @param i A character giving the origin field name in mat.
#' @param j A character giving the destination field name in mat.
#' @param fij A character giving the flow field name in mat.
#' @return A square matrix of flows. Diagonal can be filled or empty depending on data used.
#' @examples
#' # Import data
#' nav <- read.csv(system.file("csv/nav.csv", package = "flows"))
#' # Prepare data
#' myflows <- prepare_mat(x = nav, i = "i", j = "j", fij = "fij")
#' myflows[1:5,1:5]
#' @export
prepare_mat <- function(x, i, j, fij){
  mat <- x[,c(i,j,fij)]
  names(mat) <- c("i", "j", "fij")
  listUnits <- unique(c(unique(mat$i),unique(mat$j)))
  matfinal <- matrix(nrow = length(listUnits), ncol = length(listUnits),
                     dimnames = list(listUnits, listUnits))
  dmat <- reshape(mat, direction = "wide", idvar = "i", timevar = "j", sep = "_x_")
  row.names(dmat) <- dmat[, 1]
  dmat <- dmat[, -1]
  dmat <- as.matrix(dmat)
  colnames(dmat) <- unlist(lapply(
    (strsplit(colnames(dmat), split = "_x_", fixed = TRUE)),
    function(x)x[2])
  )
  i <- factor(row.names(dmat), levels = row.names(dmat))
  j <- factor(colnames(dmat), levels = colnames(dmat))
  matfinal[levels(i), levels(j)] <- dmat
  matfinal[is.na(matfinal)] <- 0
  return(matfinal)
}


