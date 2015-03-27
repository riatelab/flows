#' @title flow
#' @name flow
#' @description This package contains various functions.
#' @docType package
NULL

#' @title MRE44
#' @name MRE44
#' @references
#' http://www.insee.fr/fr/themes/detail.asp?reg_id=99&ref_id=migration-residentielle-08
#' @docType data
NULL

#' @title COM44
#' @name COM44
#' @references
#' http://professionnels.ign.fr/geofla#tab-3
#' @docType data
NULL




#' @title Flows Preparation
#' @name flowPrep
#' @description From a long format matrix to a weight data.frame and a wide
#' squared matrix
#' @param mat A data.frame of flows between \code{i} and \code{j} (long format
#' matrix: \code{(i, j, fij)}.
#' @param i A character giving the origin field in \code{mat}.
#' @param j A character giving the destination field in \code{mat}.
#' @param fij A character giving the flow field in \code{mat}.
#' @examples
#'data(LoireAtlantique)
#' @export
flowPrep <- function(mat, i, j, fij){
  mat <- mat[,c(i,j,fij)]
  names(mat) <- c("i", "j", "fij")
  listUnits <- unique(c(unique(mat$i),unique(mat$j)))
  fullMat <- expand.grid(listUnits,listUnits, stringsAsFactors = F)
  names(fullMat) <- c("i","j")
  fullMat <- merge(fullMat,mat,by=c("i","j"),all.x=TRUE)
  fullMat <- reshape2::dcast(data = fullMat, formula = i~j, value.var="fij",
                             fill = 0, sum)
  row.names(fullMat) <- fullMat[,1]
  fullMat <- fullMat[, -1]
  fullMat <- as.matrix(fullMat)
  fullMat[is.na(fullMat)] <- 0
  w <- data.frame(id = row.names(fullMat),
                  sumOut = rowSums(fullMat),
                  sumIn = colSums(fullMat))
  return(list(w = w, mat = fullMat))
}











#' @title Dominant Flows Analysis
#' @name flowDom
#' @description Computes dominant flows on a matrix :
#'
#' \code{i} is dominated by \code{j} when \code{i} sends its biggest flow to
#' \code{j} and \code{j} receives more total inbound flows than \code{i}. A
#' variant is that \code{i} is dominated by \code{j} when \code{i} sends at
#' least a specific share (\code{k}) of its outbound flows to \code{j} and
#' \code{j} receives more total inbound flows than \code{i}.
#' @param mat A data.frame of flows between from flowPrep.
#' @param w A vector of weight
#' @param k A numeric giving the minimum share of \code{i} outbound flows that
#' must be sent to \code{j}. (variant - optional)
#' @details The output of the function is a data.frame giving couples where
#' \code{i} is dominated by \code{j}. Other variables are :\itemize{
#'  \item{dom: }{share of \code{i} outbound flows sent to \code{j}.}
#'  \item{fij: }{flow from \code{i} to \code{j}}
#'  \item{sumInI: }{\code{i} total inbound flows}
#'  \item{sumOutI: }{\code{i} total outbound flows}
#'  \item{sumInJ: }{\code{j} total inbound flows}
#'  \item{sumOutJ: }{\code{j} total outbound flows}
#' }
#' @examples
#'data(LoireAtlantique)
#'dom1<- flowDom(mat = MRE44,
#'               i = "DCRAN",
#'               j = "CODGEO",
#'               fij = "NBFLUX_C08_POP05P")
#'## or
#'dom2 <- flowDom(mat = MRE44,
#'                i = "DCRAN",
#'                j = "CODGEO",
#'                fij = "NBFLUX_C08_POP05P",
#'                k = 0.5)
#' @references Nystuen, J. D. and Dacey, M. F. (1961), A GRAPH THEORY
#' INTERPRETATION OF NODAL REGIONS. Papers in Regional Science, 7: 29-42.
#' @export
flowDom <- function(mat, w, k = NULL){
#   mat <- mat[,c(i,j,fij)]
#   names(mat) <- c("i", "j", "fij")
#   listUnits <- unique(c(unique(mat$i),unique(mat$j)))
#   fullMat <- expand.grid(listUnits,listUnits, stringsAsFactors = F)
#   names(fullMat) <- c("i","j")
#   fullMat <- merge(fullMat,mat,by=c("i","j"),all.x=TRUE)
#   fullMat <- reshape2::dcast(data = fullMat, formula = i~j, value.var="fij",
#                              fill=0, sum)
#   row.names(fullMat) <- fullMat[,1]
#   fullMat <- fullMat[, -1]
#   fullMat <- as.matrix(fullMat)
#   fullMat[is.na(fullMat)] <- 0
#   diag(fullMat) <- 0
#   dimMat <- dim(fullMat)[1]
#   sumIn <- apply(fullMat, 2, sum)
#   sumOut <- apply(fullMat, 1, sum)
#   maxOut <- apply(fullMat, 1, max)
  DOM <- matrix(data = 0, nrow = dim(mat)[1], ncol = dim(mat)[1],
                dimnames = list(rownames(mat),colnames(mat) ))
  pb <- txtProgressBar(min=0, max=dimMat, style=3)
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
      setTxtProgressBar(pb, i)
    }
  } else {
    for (i in 1:dimMat){
      maxOut <- max(mat[i,])
      sumi <- w[i]
      for (j in 1:dimMat)
      {
        if (mat[i,j] == maxOut &&
            sumi < w[j] &&
            mat[i,j] != 0) {
          DOM[i,j] <- mat[i,j]/sumi
        }
      }
      setTxtProgressBar(pb, i)
    }
  }
  X <- reshape2::melt(DOM)
  names(X)<-c("i","j","dom")
  X$i<-as.character(X$i)
  X$j<-as.character(X$j)
  X <- X[X$dom > 0, ]
  X <- merge(X, mat, by =c("i","j"), all.x = T)
  sumInAndOut <- data.frame(id = colnames(fullMat), sumIn = sumIn, sumOut =
                              sumOut)
  X <- merge(X, sumInAndOut, by.x = "i", by.y = "id", all.x = T)
  X <- merge(X, sumInAndOut, by.x = "j", by.y = "id", all.x = T, suffixes =
               c('I','J' ))
  close(pb)
  return(X)
}




#' @title Plot Dominant Flows
#' @name plotflowDom
#' @description Plot a map of the Dominant Flows. It uses, as input, the output
#' of the \code{\link{flowDom}} function.
#' @param fdom A data.frame outputed by the \code{\link{flowDom}} function.
#' @param spdf A SpatialPolygonsDataFrame to be linked to \code{fdom}.
#' @param id A character giving the identifier field in \code{spdf} to be linked
#'  to \code{i} and \code{j}.
#' @param name A character giving the label field in \code{spdf} to be
#' plotted. A click on a map unit will prompt a unit name. (interactive session
#' only - optional)
#' @details The output of the function is a plot of a map. The map shows which
#' units are either "dominant", "dominated" or both. The darker a link is the
#' higher the share of the flow is in its total outbound flows.
#' @examples
#'data(LoireAtlantique)
#'dom1<- flowDom(mat = MRE44,
#'               i = "DCRAN",
#'               j = "CODGEO",
#'               fij = "NBFLUX_C08_POP05P")
#'if (interactive()){
#'  plotflowDom(fdom = dom1, spdf = COM44, id = "INSEE_COM", name = "NOM_COM")
#'}
#' @export
plotflowDom <- function(fdom, spdf, id, name = NULL ){

  sumi <- unique(fdom[,c("i","sumInI")])
  sumj <- unique(fdom[,c("j","sumInJ")])
  names(sumj) <- c("i", "sumInI")
  sumij <- unique(rbind(sumi,sumj))

  pts <- data.frame(sp::coordinates(spdf),spdf@data)
  names(pts)[1:2] <- c("X", "Y")
  pts <- merge(pts, sumij[,c("i", "sumInI")], by.x = id, by.y = "i", all.x = T)

  fdom <- merge(fdom, pts, by.x = "i", by.y = id, all.x = T,
                suffixes = c("i","j"))
  fdom <- merge(fdom, pts, by.x = "j", by.y = id, all.x = T,
                suffixes = c("i","j"))
  listI <- unique(fdom$i)
  listJ <- unique(fdom$j)

  x <- data.frame(pts, d = pts[,id] %in% listI, D = pts[,id] %in% listJ,
                  col = NA, pch = NA, cex = pts$sumInI)

  bbbox <- sp::bbox(spdf)
  x1 <- bbbox[1]
  y1 <- bbbox[2]
  x2 <- bbbox[3]
  y2 <- bbbox[4]
  sfdc <- (x2-x1)*(y2-y1)
  sc <- sum(x$sumInI,na.rm=TRUE)
  x$cex <- sqrt((x$sumInI * 0.02 * sfdc / sc) / pi)
  x <- x[order(x$cex,decreasing=TRUE),]
  x[!x$d & x$D, "col"] <- "#cc2a36"
  x[x$d & x$D, "col"] <- "#eb6841"
  x[x$d & !x$D,"col"] <- "#edc951"

  cl <- seq(min(fdom$dom), max(fdom$dom),length.out = 4)
  fdom$col <- findInterval(fdom$dom, cl,all.inside = T)

  fdom[fdom$col == 1,"col"] <- "#4d343440"
  fdom[fdom$col == 2,"col"] <- "#4d343470"
  fdom[fdom$col == 3,"col"] <- "#4d3434"
  #   fdom[fdom$dom <= .1,"col"] <- "#4d343410"
  #   fdom[fdom$dom > .1,"col"] <-
  #     paste("#4d3434", round(fdom[fdom$dom >.1,"dom"]*100), sep="")
  #   fdom[fdom$dom >= 1,"col"] <- "#4d3434"

  fdom <- fdom[order(-fdom$dom),]

  sp::plot(spdf, col = "#cceae7", border = "grey70")

  segments(fdom$Xi, fdom$Yi, fdom$Xj, fdom$Yj, col=fdom$col, lwd = 4)

  symbols(x[,c("X","Y")], circles = x$cex, add = TRUE, bg = x$col,
          fg ="grey50",
          inches = FALSE)

  legend(x = "topleft",
         legend =  c("Dominant", "Dominant-Dominated", "Dominated",
                     "Share of the sent flows :",
                     paste("high (",round(cl[3]*100,0),"-",round(cl[4]*100,0),
                           ")",sep = ""),
                     paste("medium (",round(cl[2]*100,0),"-",round(cl[3]*100,0),
                           ")",sep = ""),
                     paste("low (",round(cl[1]*100,0),"-",round(cl[2]*100,0),
                           ")",sep = "")),
         col= c(rep("grey50",3),rep(NA,4)),
         cex = 0.7,
         pch = c(21,21,21,16,22,22,22),
         pt.bg = c("#cc2a36", "#eb6841","#edc951","#ffffff00","#4d3434",
                   "#4d343470","#4d343440"),
         pt.cex = c(2,1,0.5,2,2,2,2))
  if (!is.null(name)){
    if(interactive()){
      x <- locator()
      if (!is.null(x)){
        X <- sp::SpatialPoints(data.frame(x), proj4string = spdf@proj4string)
        text(x = x$x, y = x$y, labels = sp::over(X, spdf)[,name], cex = 0.6,
             adj = c(0,0))
      }
    } else {
      text(x[x$col == "#cc2a36",c("X","Y")],labels = x[x$col == "#cc2a36", name],
           cex = 0.6)
    }
  }
}




