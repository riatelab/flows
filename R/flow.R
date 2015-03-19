#' @title flow
#' @name flow
#' @description This package contains various functions.
#' @docType package
NULL

#' @title MRE46
#' @name MRE46
#' @references
#' http://www.insee.fr/fr/themes/detail.asp?reg_id=99&ref_id=migration-residentielle-08
#' @docType data
NULL

#' @title COM46
#' @name COM46
#' @references
#' http://professionnels.ign.fr/geofla#tab-3
#' @docType data
NULL

#' @title Compute dominant flows
#' @name flowDom
#' @description Compute dominant flows on a long matrix (i, j, fij). Two
#' methods are allowed: i sends its major flow to j or i sends a specific share
#' to j.
#' @param mat A data.frame of flows between \code{i} and \code{j} (long format).
#' @param i A character giving the name of emiting units field in \code{mat}.
#' @param j A character giving the name of receiving units field in \code{mat}.
#' @param fij A character giving the name of the flow field in \code{mat}.
#' @param k A numeric giving the minimum share of flow that must be emitted from
#'  \code{i} to \code{j}. (optional)
#' @details The output of the function is a data.frame giving couples where
#' \code{i} is dominated by \code{j} and the share of outputs from i to j in the
#' total outputs from i.
#' @examples
#'data(lot)
#'dom1<- flowDom(mat = MRE46,
#'               i = "DCRAN",
#'               j = "CODGEO",
#'               fij = "NBFLUX_C08_POP05P")
#'dom2 <- flowDom(mat = MRE46,
#'                i = "DCRAN",
#'                j = "CODGEO",
#'                fij = "NBFLUX_C08_POP05P",
#'                k = 0.5)
#' @references Nystuen, J. D. and Dacey, M. F. (1961), A GRAPH THEORY
#' INTERPRETATION OF NODAL REGIONS. Papers in Regional Science, 7: 29-42.
#' @export
flowDom <- function(mat, i, j, fij, k = NULL){
  mat <- mat[,c(i,j,fij)]
  names(mat) <- c("i", "j", "fij")
  listUnits <- unique(c(unique(mat$i),unique(mat$j)))
  fullMat <- expand.grid(listUnits,listUnits, stringsAsFactors = F)
  names(fullMat) <- c("i","j")
  fullMat <- merge(fullMat,mat,by=c("i","j"),all.x=TRUE)
  fullMat <- reshape2::dcast(data = fullMat, formula = i~j, value.var="fij",
                             fill=0, sum)
  row.names(fullMat) <- fullMat[,1]
  fullMat <- fullMat[, -1]
  fullMat <- as.matrix(fullMat)
  fullMat[is.na(fullMat)] <- 0
  diag(fullMat) <- 0
  dimMat <- dim(fullMat)[1]
  sumIn <- apply(fullMat, 2, sum)
  sumOut <- apply(fullMat, 1, sum)
  maxOut <- apply(fullMat, 1, max)
  DOM <- matrix(data = 0, nrow = dimMat, ncol = dimMat,
                dimnames = list(rownames(fullMat),colnames(fullMat) ))
  if(!is.null(k)){
    for (i in 1:dimMat){
      for (j in 1:dimMat)
      {
        if ((fullMat[i,j] / sumOut[i]) >= k &&
            sumIn[i] < sumIn[j] &&
            fullMat[i,j] != 0) {
          DOM[i,j] <- (fullMat[i,j] / sumOut[i])
        }
      }
    }
  } else {
    for (i in 1:dimMat){
      for (j in 1:dimMat)
      {
        if (fullMat[i,j] == maxOut[i] &&
            sumIn[i] < sumIn[j] &&
            fullMat[i,j] != 0) {
          DOM[i,j] <- (fullMat[i,j] / sumOut[i])
        }
      }
    }
  }
  X <- reshape2::melt(DOM)
  names(X)<-c("i","j","dom")
  X$i<-as.character(X$i)
  X$j<-as.character(X$j)
  X <- X[X$dom > 0, ]
  return(X)
}


#' @title Plot dominant flows
#' @name plotflowDom
#' @description Plot a map of the Dominant Flows. It uses, as input, the output
#' of the flowDom function.
#' @param fdom A data.frame outputed by the flowDom function.
#' @param spdf A SpatialPolygonsDataFrame to be linked to fdom.
#' @param id A character giving the name of the 'id' field in spdf.
#' @param name A character giving the name of a label field in spdf to be
#' plotted. (optional)
#' @details The output of the function is a plot of a map. The map shows which
#' units are either "dominant", "dominated" or both. The darker a link is the
#' higher the share of the flow is in the total flows emitted by the unit.
#' @examples
#'data(lot)
#'dom1<- flowDom(mat = MRE46,
#'               i = "DCRAN",
#'               j = "CODGEO",
#'               fij = "NBFLUX_C08_POP05P")
#'dom2 <- flowDom(mat = MRE46,
#'                i = "DCRAN",
#'                j = "CODGEO",
#'                fij = "NBFLUX_C08_POP05P",
#'                k = 0.5)
#'plotflowDom(fdom = dom1, spdf = COM46, id = "INSEE_COM", name = "NOM_COM")
#' @export
plotflowDom <- function(fdom, spdf, id, name = NULL ){
  pts <- data.frame(sp::coordinates(spdf),spdf@data)
  names(pts)[1:2] <- c("X", "Y")
  fdom <- merge(fdom, pts, by.x = "i", by.y = id, all.x = T,
                suffixes = c("i","j"))
  fdom <- merge(fdom, pts, by.x = "j", by.y = id, all.x = T,
                suffixes = c("i","j"))

  listI <- unique(fdom$i)
  listJ <- unique(fdom$j)
  x <- data.frame(pts, d = pts$INSEE_COM %in% listI, D = pts$INSEE_COM %in% listJ,
                  col = NA, pch = NA, cex = NA)
  # D
  x[!x$d & x$D, "col"] <- "#cc2a36"
  x[!x$d & x$D, "cex"] <- 2
  # D d
  x[x$d & x$D, "col"] <- "#eb6841"
  x[x$d & x$D, "cex"] <- 1
  # d
  x[x$d & !x$D,"col"] <- "#edc951"
  x[x$d & !x$D, "cex"] <- 0.5

  cl <- seq(min(fdom$dom), max(fdom$dom),length.out = 4)
  fdom$col <- findInterval(fdom$dom, cl,all.inside = T)
  fdom[fdom$dom != 1,"col"] <- paste("#4d3434", round(fdom[fdom$dom!=1,"dom"]*100),sep="")
  fdom[fdom$dom == 1,"col"] <- "#4d3434"
  # fdom[,"col"] <- paste("#4d3434", round(fdom$dom*100)}else{99},sep="")
  fdom <- fdom[order(fdom$dom),]

  sp::plot(spdf, col = "#cceae7", border = "grey70")
  segments(fdom$Xi, fdom$Yi, fdom$Xj, fdom$Yj, col=fdom$col, lwd = 4)

  points(x[,1:2], pch = 21, col = "grey50", bg = x$col, cex = x$cex)
  legend(x = "topleft",
         legend =  c("Dominant", "Dominant-Dominated", "Dominated",
                     "Share of the emitted flow :",
                     paste("higher - max = ",round(max(cl),2),sep = ""),
                     paste("lower - min = ", round(min(cl),2),sep = "")),
         col= c(rep("grey50",3),rep(NA,3)),
         cex = 0.7,
         pch = c(21,21,21,16,22,22),
         pt.bg = c("#cc2a36", "#eb6841","#edc951","#ffffff00","#4d3434",
                   "#4d343420"),
         pt.cex = c(2,1,0.5,0.1,2,2))
  if (!is.null(name)){
    text(x[x$cex == 2,1:2],labels = x[x$cex == 2, name], cex = 0.6)
  }
}




