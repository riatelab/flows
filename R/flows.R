#' @title flows
#' @name flows
#' @description This package contains various functions to filter flows matrices (dominant and major flows), provides statistics on selection made
#'  and proposes map and graph visualisation.
#' @docType package
NULL

#' @title MRE44
#' @name MRE44
#' @description  Migration flows at the communal level for departement of Loire Atlantique (change of residence between 2008 and 2012)
#' @details
#' Variables: \cr
#' \itemize{
#' \item{CODGEO: Departement and commune of residence}
#' \item{LIBGEO: Name of the commune of residence}
#' \item{DCRAN: INSEE code of former departement and commune of residence (municipal arrondissement for Paris, Lyon, Marseille)}
#' \item{L_DCRAN: Name of former commune (municipal arrondissement for Paris, Lyon, Marseille)}
#' \item{NBFLUX_C08_POP05P: Number of migrants (5 years or more) from commune or arrondissement DCRAN to commune CODGEO}}
#' @references
#' \url{http://www.insee.fr/fr/themes/detail.asp?reg_id=99&ref_id=migration-residentielle-08}
#' @docType data
NULL

#' @title COM44
#' @name COM44
#' @description SpatialPolygonsDataFrame of Loire Atlantique (communal level) - same variables than \link{MRE44}
#' @references
#' \url{http://professionnels.ign.fr/geofla#tab-3}
#' @docType data
NULL


#' @title MRE31
#' @name MRE31
#' @description Migration flows at the communal level for departement of Haute Garonne (change of residence between 2008 and 2012)
#' @details
#' Variables: \cr
#' \itemize{
#' \item{CODGEO: Departement and commune of residence}
#' \item{LIBGEO: Name of the commune of residence}
#' \item{DCRAN: INSEE code of former departement and commune of residence (municipal arrondissement for Paris, Lyon, Marseille)}
#' \item{L_DCRAN: Name of former commune (municipal arrondissement for Paris, Lyon, Marseille)}
#' \item{NBFLUX_C08_POP05P: Number of migrants (5 years or more) from commune or arrondissement DCRAN to commune CODGEO}}
#' @references
#' \url{http://www.insee.fr/fr/themes/detail.asp?reg_id=99&ref_id=migration-residentielle-08}
#' @docType data
NULL

#' @title COM31
#' @name COM31
#' @description SpatialPolygonsDataFrame of Haute Garonne (communal level) - same variables than \link{MRE31}
#' @references
#' \url{http://professionnels.ign.fr/geofla#tab-3}
#' @docType data
NULL

#### Public

#' @title Flows Preparation
#' @name prepflows
#' @description From a long format matrix to a data.frame with sum of outputs and sum of inputs and a wide
#' matrix of flows.
#' @param mat A data.frame of flows between origins and destinations: long format
#' matrix (origins, destinations, flows intensity).
#' @param i A character giving the origin field name in \code{mat}.
#' @param j A character giving the destination field name in \code{mat}.
#' @param fij A character giving the flow field name in \code{mat}.
#' @return A square flows matrix. Diagonal can be filled or empty depending on data used.
#' @examples
#' data(LoireAtlantique)
#' myflows <- prepflows(mat = MRE44, i = "DCRAN", j = "CODGEO", fij = "NBFLUX_C08_POP05P")
#' myflows[1:5,1:5]
#' @import reshape2
#' @export
prepflows <- function(mat, i, j, fij){
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
  #   w <- data.frame(id = row.names(fullMat),
  #                   sumOut = rowSums(fullMat),
  #                   sumIn = colSums(fullMat))
  # return(list(dfw = w, mat = fullMat))
  return(fullMat)
}

#' @title Descriptive Statistics of Flows Matrix
#' @name statmat
#' @description Give various indicators and graphical outputs on flow matrices
#' @param mat A square flow matrix.
#' @return  The function returns graphics, statistics and a list. \cr
#' \itemize{
#' \item{nblinks: number of cells with values > 0,}
#' \item{density: number of links / number of possible links (also called gamma index by geographers)}
#' \item{connectcomp: number of connected components (isolates included,
#' weakly connected see \code{\link{clusters}})}
#' \item{connectcompx: number of connected components (isolates deleted,
#' weakly connected see \code{\link{clusters}})}
#' \item{sizecomp: a data.frame only returned in the list, size of the connected components.}
#' \item{compocomp: a data.frame only returned in the list, connected components membership of units.}
#' \item{sumflows: sum of flows}
#' \item{min: min flow }
#' \item{Q1: Q1 flow}
#' \item{median: median flow}
#' \item{Q3: Q3 flow}
#' \item{max: max flow}
#' \item{mean: mean flow}
#' \item{sd: standart deviation flow}}
#' Graphics returned are: \cr
#' \itemize{
#' \item{degree distribution}
#' \item{weigthed degree distribution}
#' \item{Lorenz curve regarding links weights}
#' \item{Boxplot regarding links weights}
#' }
#' @import igraph
#' @examples
#' data(LoireAtlantique)
#'myflows <- prepflows(mat = MRE44, i = "DCRAN", j = "CODGEO", fij = "NBFLUX_C08_POP05P")
#'x <- statmat(myflows)
#' # Size of connected components
#' x$sizecomp
#' # Sum of flows
#' x$sumflows
#' x <- statmat(myflows)
#' # Size of connected components
#' x$sizecomp
#' # Sum of flows
#' x$sumflows
#' @export
statmat <- function(mat){
  nbcell <- length(mat)
  matdim <- dim(mat)[1]
  sumflows <- sum(mat)
  matbool <- mat
  matbool[mat > 0] <- 1
  nbcellfull <- sum(matbool)
  vmat <- as.vector(mat[mat > 0])
  vmat <- vmat[order(vmat, decreasing = FALSE)]
  sumflows <- sum(vmat)
  vmatcs <- cumsum(vmat) /  sumflows * 100
  summaryflows <- summary(vmat)
  summaryflows <- c(summaryflows, sd(vmat))
  names(summaryflows) <- NULL

  ## graphic outputs
  old.par <- par (mfrow = c(2,2))
  ## rank-size link
  deg <- rowSums(matbool)
  deg <- deg[deg>0]
  plot(deg[order(deg, decreasing = TRUE)], type = "l", log = "xy",
       xlab = "rank (log)", ylab = "size (log nb. flows)")
  title("rank - size")
  ## rank size flow
  deg <- rowSums(mat)
  deg <- deg[deg>0]
  plot(deg[order(deg, decreasing = TRUE)], type = "l", log = "xy",
       xlab = "rank (log)", ylab = "size (log flow intensity)")
  title("rank - size (weighted)")
  ##lorenz
  plot( y = vmatcs, x = seq(0,100,length.out = length(vmatcs)), type = "l",
        xlim = c(0,100), ylim = c(0,100),
        xlab = "cum. nb. flows", ylab = "cum. intensity of flows")
  title ("Lorenz Curve")
  ## boxplot
  boxplot(as.vector(mat[mat>0]), log = "y")
  title("Boxplot")
  par(old.par)

  ## Connected components of a graph
  g <- graph.adjacency(adjmatrix = mat, mode = "directed", weighted = TRUE)
  clustg <- clusters(graph = g, mode = "weak")

  connectcomp <- clustg$no
  connectcompx <- length(clustg$csize[clustg$csize>1])
  sizecomp <- data.frame(idcomp = seq(1, length(clustg$csize)),
                          sizecomp = clustg$csize)
  compocomp <-  data.frame(id = V(g)$name, idcomp = clustg$membership)

  ## stat cat
  cat('matrix dimension:', matdim, "X", matdim,"\n" )
  cat('nb. links:', nbcellfull, "\n" )
  cat('density:', nbcellfull/nbcell, "\n" )
  cat('nb. of components (weak)', connectcomp, "\n")
  cat("nb. of components (weak, size > 1)", connectcompx, "\n")
  cat('sum of flows:', sumflows, "\n")
  cat('min:', summaryflows[1] ,"\n")
  cat('Q1:', summaryflows[2] ,"\n")
  cat('median:', summaryflows[3] ,"\n")
  cat('Q3:', summaryflows[5] ,"\n")
  cat('max:', summaryflows[6] ,"\n")
  cat('mean:', summaryflows[4] ,"\n")
  cat('sd:', summaryflows[7] ,"\n")

  ## stat list
  matstat <- list(matdim = dim(mat),
                  nblinks = nbcellfull,
                  density = nbcellfull/nbcell,
                  connectcomp = connectcomp,
                  connectcompx = connectcompx,
                  sizecomp = sizecomp,
                  compocomp = compocomp,
                  sumflows = sumflows,
                  min =  summaryflows[1],
                  Q1 =  summaryflows[2],
                  median = summaryflows[3],
                  Q3 = summaryflows[5],
                  max = summaryflows[6],
                  mean = summaryflows[4],
                  sd = summaryflows[7]
  )
  return(invisible(matstat))

}



#' @title Flow Selection From Origins
#' @name firstflows
#' @description Flow selection from \emph{i}.
#' @param mat A square matrix of flows \emph{fij}
#' @param method One of "nfirst", "xfirst" or "xsumfirst". \cr nfirst = select \emph{k} first \emph{fij}
#' from \emph{i}
#' \cr xfirst = select x \emph{fij} from i where \emph{fij} > k  \cr xsumfirst = select x \emph{fij} from \emph{i}
#'  while sum(\emph{fij}) < k.
#' @param ties.method In case of equality with 'nfirst' method (see \link{rank},
#' "random" and "first" select only one flow, other methods may select various).
#' @param k Selection threshold (can be relative or absolute).
#' @references J. Nystuen & M. Dacey, 1961, A graph theory interpretation of nodal flows,
#' \emph{Papers and Proceedings of the Regional Science Association}, vol. 7,  29-42.
#' @return A boolean matrix of selected flows. To incorporate flows intensity, making the product
#' mat * boolean matrix is necessary.
#' @examples
#' data(LoireAtlantique)
#' myflows <- prepflows(mat = MRE44, i = "DCRAN", j = "CODGEO", fij = "NBFLUX_C08_POP05P")
#'
#' #remove diagonal
#' diag(myflows) <- 0
#' statmat(myflows)
#'
#' #select 2 flows per spatial unit
#' fflows1 <- firstflows(myflows, method = "nfirst", ties.method = "first", 2)
#' fflow1 <- fflows1 * myflows
#' statmat(fflow1)
#'
#' #select flows > 20
#' fflows2 <- firstflows(myflows, method = "xfirst", ties.method = "first", 20)
#' fflow2 <- fflows2 * myflows
#' statmat(fflow2)
#'
#' #select sum(flows) from i > 20
#' fflows3 <- firstflows(myflows, method = "xsumfirst", ties.method = "first", 20)
#' fflow3 <- fflows3 * myflows
#' statmat(fflow3)
#'
#' # Select each flows that represent at least 10% of the outputs
#' myflowspct <- myflows / rowSums(myflows) * 100
#' fflows2 <- firstflows(mat = myflowspct, method = "xfirst", ties.method = "first", k = 10)
#' fflow2 <- fflows2 * myflows
#' statmat(fflow2)
#' @export
firstflows <- function(mat, method = "nfirst", ties.method = "first",k){
  # list of i, j selected
  lfirst <- apply(mat, 1, get(method), k = k, ties.method = ties.method)
  # if only one selected
  if(is.null(dim(lfirst))){
    lfirst <- as.list((lfirst))
  }
  # control class output
  if(is.matrix(lfirst)) {
    lfirst <- as.list(as.data.frame(lfirst, stringsAsFactors = FALSE))
  }
  matfinal <- mat
  matfinal[] <- 0
  # control 0 selection
  if (length(lfirst)<1){
    return(matfinal)
  }
  for (i in 1:nrow(matfinal)){
    matfinal[names(lfirst[i][]),lfirst[[i]]] <- 1
  }
  return(matfinal)
}


#' @title Flow Selection From the Total Matrix
#' @name firstflowsg
#' @description Various flow selection on global criterions
#' @param mat A square matrix of flows \emph{fij}
#' @param method One of "nfirst", "xfirst" or "xsumfirst". \cr nfirst = select \emph{k} first \emph{fij}
#' from the input matrix
#' \cr xfirst = select x \emph{fij} from the input matrix where \emph{fij} > k  \cr xsumfirst = select x \emph{fij} from the input matrix
#'  while sum(\emph{fij}) < k.
#' @param k Selection threshold (can be relative or absolute).
#' @param ties.method In case of equality with 'nfirst' method (see \link{rank},
#' "random" and "first" select only one flow, other methods may select various).
#' @return A boolean matrix of selected flows. To incorporate flows intensity, making the product
#' mat * boolean matrix is necessary.
#' @examples
#' data(LoireAtlantique)
#' myflows <- prepflows(mat = MRE44, i = "DCRAN", j = "CODGEO", fij = "NBFLUX_C08_POP05P")

#' #remove diagonal
#' diag(myflows) <- 0
#' statmat(myflows)

#' #select 50 first flow on the whole matric

#' fflows1 <- firstflowsg(myflows, method = "nfirst", ties.method = "first", 50)
#' fflow1 <- fflows1 * myflows
#' statmat(fflow1)

#' #select flows > 50

#' fflows2 <- firstflowsg(myflows, method = "xfirst", ties.method = "first", 50)
#' fflow2 <- fflows2 * myflows
#' statmat(fflow2)

#' #select sum(flows) > 50% of the whole matrix
#' #I don't know how to make a relative selection...
#' @export
firstflowsg <- function(mat, method = "nfirst", k, ties.method = "first"){
  matfinal <- mat
  matfinal[] <- 0
  if (method == "nfirst"){
    matfinal[rank(mat, ties.method = ties.method) > ((dim(mat)[1]*dim(mat)[2]) - k)] <- 1
  }
  if (method == "xfirst"){
    matfinal[mat >= k] <- 1
  }
  if (method == "xsumfirst"){
    matv <- as.vector(mat)
    names(matv) <- 1:length(matv)
    matvo <- matv[order(matv, decreasing = TRUE)]
    matvo <- cumsum(matvo)
    nbgood <- (length(matvo[matvo < k ])+1)
    matvo[] <- c(rep(1,nbgood), rep(0,(length(matvo)-nbgood)))
    matfinal[] <- matvo[order(as.numeric(names(matvo)), decreasing = FALSE)]
  }
  matfinal[mat == 0] <- 0
  return(matfinal)
}


#' @title Dominant Flows Selection
#' @name domflows
#' @description Compute the dominant flow analysis based on the paper from Nystuen and Dacey
#' @param mat A square matrix of flows
#' @param wi A vector of weight for i
#' @param wj A vector of weight for j
#' @param k Threshold. wj/wi> k
#' @references J. Nystuen & M. Dacey, 1961, A graph theory interpretation of nodal flows,
#' \emph{Papers and Proceedings of the Regional Science Association}, vol. 7,  29-42.
#' @return A boolean matrix of selected flows. To incorporate flows intensity, making the product
#' mat * boolean matrix is necessary.
#' @examples
#' data(LoireAtlantique)
#' myflows <- prepflows(mat = MRE44, i = "DCRAN", j = "CODGEO", fij = "NBFLUX_C08_POP05P")

#' #remove diagonal
#' diag(myflows) <- 0
#' statmat(myflows)

#' #create the tree of dominant, intermediary and dominated nodes
#' domf <- domflows(mat = myflows, wi = colSums(myflows), wj = colSums(myflows), k = 1)
#' domfw <-domf * myflows
#' statmat(domfw)
#' @export
domflows <- function(mat, wi, wj, k){
  # list of i, j selected
  matfinal <- mat
  matfinal[] <- 0
  for (i in 1:dim(mat)[1]){
    for (j in 1:dim(mat)[2]){
      if (wi[i] > 0){
        if ((wj[j]/wi[i]) > k){
          matfinal[i,j] <- 1
        }
      }
    }
  }
  matfinal[mat == 0] <- 0
  return(matfinal)
}

#' @title Dominant Flows Graph
#' @name plotDomFlows
#' @description Display a dominant flows graph
#' @param mat A square matrix of flows
#' @note We do not propose visualisation for other outputs as square matrices produced can easily be plot
#' (after transformation) with \link[igraph]{plot.igraph} or \link[sna]{gplot} functions.
#' @examples
#' data(LoireAtlantique)
#' mat <- prepflows(mat = MRE44, i = "DCRAN", j = "CODGEO", fij = "NBFLUX_C08_POP05P")
#' diag(mat) <- 0
#' x <- domflows(mat = mat, wi = colSums(mat), wj = colSums(mat), k = 1)
#' firstx <- firstflows(mat = mat, method = "nfirst", ties.method = "first", k = 1)
#' xnb <- firstflows(mat = mat, method = "xfirst", ties.method = "first", k = 20)
#' mat <- mat * firstx * x * xnb
#' plotDomFlows(mat)
#' @export
plotDomFlows <- function(mat){
  g <- graph.adjacency(adjmatrix = mat,mode = "directed", weighted = TRUE)
  g <- delete.vertices(g, names(degree(g)[degree(g)==0]))
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
#   lg <- layout.fruchterman.reingold(g)
#   g <- set.graph.attribute(graph = g, name = "layout", value = lg)

  x <- igraph::plot.igraph(g, vertex.label = V(g)$names, vertex.label.cex = 1,
                   vertex.label.color = "black",
                   vertex.size = V(g)$size, edge.arrow.size = 0)
  # legend(x = "bottomleft", legend = c("héhé", "hoho", "haha"), col = c("red", "orange", "yellow"), pch = "-", pt.cex = 15)
  title("Dominant Flows Graph")
}

#' @title Dominant Flows Map
#' @name plotMapDomFlows
#' @description Display a dominant flows map
#' @param mat A square matrix of dominant flows
#' @param spdf Spatial*DataFrame of units
#' @param spdfid Name of the 'id' variable in the spdf data.frame
#' @param w A data.frame which contains a weight variable for the units
#' @param wid Identifier column  in w
#' @param wvar Weight variable
#' @param wcex Size factor for circles of weight
#' @param add Whether to add the plot to an existing one or not.
#' @import sp
#' @examples
#' data(LoireAtlantique)
#' mat <- prepflows(mat = MRE44, i = "DCRAN", j = "CODGEO", fij = "NBFLUX_C08_POP05P")
#' diag(mat) <- 0
#' x <- domflows(mat = mat, wi = colSums(mat), wj = colSums(mat), k = 1)
#' firstx <- firstflows(mat = mat, method = "nfirst", ties.method = "first", k = 1)
#' hab <- mat * firstx * x
#' inflows <- data.frame(id = colnames(mat), w = colSums(mat))
#' sp::plot(COM44, col = "#cceae7")
#' plotMapDomFlows(mat = hab,
#'                 spdf = COM44,
#'                 spdfid = "INSEE_COM",
#'                 w = inflows,
#'                 wid = "id",
#'                 wvar = "w",
#'                 wcex = 0.05,
#'                 add = TRUE)
#' @export
plotMapDomFlows <- function(mat, spdf,
                            spdfid, w,
                            wid, wvar,
                            wcex = 0.05, add = FALSE){
  # points management
  pts <- data.frame(sp::coordinates(spdf), id  = spdf@data[,spdfid])
  names(pts)[1:2] <- c("long", "lat")
  w <- w[,c(wid, wvar)]
  names(w) <- c("id", "var")
  pts <- merge(pts, w, by.x = "id", by.y = "id", all.x = T)

  # points size
  bbbox <- sp::bbox(spdf)
  x1 <- bbbox[1]
  y1 <- bbbox[2]
  x2 <- bbbox[3]
  y2 <- bbbox[4]
  sfdc <- (x2-x1)*(y2-y1)
  sc <- sum(pts$var, na.rm=TRUE)
  pts$cex <- sqrt((pts$var * wcex * sfdc / sc) / pi)
  pts <- pts[order(pts$cex,decreasing=TRUE),]
  pts <- pts[pts$cex > 0, ]

  # Segment management
  colnames(mat) <- paste("X", colnames(mat), sep="")
  row.names(mat) <- paste("X", row.names(mat), sep="")
  fdom <- reshape2::melt(mat)
  names(fdom) <- c("i", "j", "fij")
  fdom <- fdom[fdom$fij > 0,]
  fdom$i <- substr(x = fdom$i, 2 , nchar(as.character(fdom$i)))
  fdom$j <- substr(x = fdom$j, 2 , nchar(as.character(fdom$j)))
  fdom <- merge(fdom, pts, by.x = "i", by.y = "id", all.x = T,
                suffixes = c("i","j"))
  fdom <- merge(fdom, pts, by.x = "j", by.y = "id", all.x = T,
                suffixes = c("i","j"))
  fdom$width <- (fdom$fij * 8 / (max(fdom$fij) - min(fdom$fij))) + 2

  # points color
  pts$col <- "red"
  pts[pts$id %in% fdom$j & pts$id %in% fdom$i, "col"] <- "orange"
  pts[!pts$id %in% fdom$j & pts$id %in% fdom$i, "col"] <- "yellow"

  # Affichage points and segments
  if(add == FALSE){
    sp::plot(spdf, col = NA, border = NA, add = F)
  }
  segments(fdom$longi, fdom$lati, fdom$longj, fdom$latj,
           col="grey20", lwd = fdom$width)
  symbols(pts[,c("long","lat")],
          circles = pts$cex,
          add = TRUE,
          bg = pts$col,
          fg ="grey50",
          inches = F)
  segments(fdom$longi, fdom$lati, fdom$longj, fdom$latj,
           col="#00000010", lwd = fdom$width)

  # Affichage legend
  legend( x="topleft",
          legend=round(c(min(fdom$fij),max(fdom$fij)),0),
          col=c("grey20"), lwd=c(2,8), lty=1)
}









#### Private
#' @title nfirst
#' @name nfirst
#' @noRd
nfirst <- function(x, k, ties.method){
  x <- x[x > 0]
  if (length(x) > 0){
    if (length(x) > k ){
      x <- x[rank(x, ties.method = ties.method) > (length(x) - k)]
      x <- names(x)
    }else{
      x <- names(x)
    }
  }
  return(x)
}

#' @title xfirst
#' @name xfirst
#' @noRd
xfirst <- function(x, k, ties.method){
  x <- x[x > 0]
  if (length(x) > 0){
    x <- x[x > k]
    x <- names(x)
  }else{
    x <- names(x)
  }
  return(x)
}

#' @title xsumfirst
#' @name xsumfirst
#' @noRd
xsumfirst <- function(x, k, ties.method){
  x <- x[x > 0]
  x <- x[order(x, decreasing = TRUE)]
  x <- cumsum(x = x)
  if (length(x) > 0){
    if (x[1] >= k){
      x <- names(x[1])
    }else{
      # at least k
      x <- x[1:(length(x[(x <= k)==TRUE]) + 1)]
      # less than k
      # x <- x[(x <= k)]
      x <- names(x)
    }
  }else{
    x <- names(x)
  }
  # pb return error if k > sum(mat) or k > sum(mat/w)
  return(x)
}







# #' @title Plot Dominant Flows
# #' @name plotflowDom
# #' @description Plot a map of the Dominant Flows. It uses, as input, the output
# #' of the \code{\link{flowDom}} function.
# #' @param fdom A data.frame outputed by the \code{\link{flowDom}} function.
# #' @param spdf A SpatialPolygonsDataFrame to be linked to \code{fdom}.
# #' @param id A character giving the identifier field in \code{spdf} to be linked
# #'  to \code{i} and \code{j}.
# #' @param name A character giving the label field in \code{spdf} to be
# #' plotted. A click on a map unit will prompt a unit name. (interactive session
# #' only - optional)
# #' @details The output of the function is a plot of a map. The map shows which
# #' units are either "dominant", "dominated" or both. The darker a link is the
# #' higher the share of the flow is in its total outbound flows.
# #' @examples
# #'data(LoireAtlantique)
# #'dom1<- flowDom(mat = MRE44,
# #'               i = "DCRAN",
# #'               j = "CODGEO",
# #'               fij = "NBFLUX_C08_POP05P")
# #'if (interactive()){
# #'  plotflowDom(fdom = dom1, spdf = COM44, id = "INSEE_COM", name = "NOM_COM")
# #'}
# #' @export
# plotflowDom <- function(fdom, spdf, id, name = NULL ){
#
#   sumi <- unique(fdom[,c("i","sumInI")])
#   sumj <- unique(fdom[,c("j","sumInJ")])
#   names(sumj) <- c("i", "sumInI")
#   sumij <- unique(rbind(sumi,sumj))
#
#   pts <- data.frame(sp::coordinates(spdf),spdf@data)
#   names(pts)[1:2] <- c("X", "Y")
#   pts <- merge(pts, sumij[,c("i", "sumInI")], by.x = id, by.y = "i", all.x = T)
#
#   fdom <- merge(fdom, pts, by.x = "i", by.y = id, all.x = T,
#                 suffixes = c("i","j"))
#   fdom <- merge(fdom, pts, by.x = "j", by.y = id, all.x = T,
#                 suffixes = c("i","j"))
#   listI <- unique(fdom$i)
#   listJ <- unique(fdom$j)
#
#   x <- data.frame(pts, d = pts[,id] %in% listI, D = pts[,id] %in% listJ,
#                   col = NA, pch = NA, cex = pts$sumInI)
#
#   bbbox <- sp::bbox(spdf)
#   x1 <- bbbox[1]
#   y1 <- bbbox[2]
#   x2 <- bbbox[3]
#   y2 <- bbbox[4]
#   sfdc <- (x2-x1)*(y2-y1)
#   sc <- sum(x$sumInI,na.rm=TRUE)
#   x$cex <- sqrt((x$sumInI * 0.02 * sfdc / sc) / pi)
#   x <- x[order(x$cex,decreasing=TRUE),]
#   x[!x$d & x$D, "col"] <- "#cc2a36"
#   x[x$d & x$D, "col"] <- "#eb6841"
#   x[x$d & !x$D,"col"] <- "#edc951"
#
#   cl <- seq(min(fdom$dom), max(fdom$dom),length.out = 4)
#   fdom$col <- findInterval(fdom$dom, cl,all.inside = T)
#
#   fdom[fdom$col == 1,"col"] <- "#4d343440"
#   fdom[fdom$col == 2,"col"] <- "#4d343470"
#   fdom[fdom$col == 3,"col"] <- "#4d3434"
#   #   fdom[fdom$dom <= .1,"col"] <- "#4d343410"
#   #   fdom[fdom$dom > .1,"col"] <-
#   #     paste("#4d3434", round(fdom[fdom$dom >.1,"dom"]*100), sep="")
#   #   fdom[fdom$dom >= 1,"col"] <- "#4d3434"
#
#   fdom <- fdom[order(-fdom$dom),]
#
#   sp::plot(spdf, col = "#cceae7", border = "grey70")
#
#   segments(fdom$Xi, fdom$Yi, fdom$Xj, fdom$Yj, col=fdom$col, lwd = 4)
#
#   symbols(x[,c("X","Y")], circles = x$cex, add = TRUE, bg = x$col,
#           fg ="grey50",
#           inches = FALSE)
#
#   legend(x = "topleft",
#          legend =  c("Dominant", "Dominant-Dominated", "Dominated",
#                      "Share of the sent flows :",
#                      paste("high (",round(cl[3]*100,0),"-",round(cl[4]*100,0),
#                            ")",sep = ""),
#                      paste("medium (",round(cl[2]*100,0),"-",round(cl[3]*100,0),
#                            ")",sep = ""),
#                      paste("low (",round(cl[1]*100,0),"-",round(cl[2]*100,0),
#                            ")",sep = "")),
#          col= c(rep("grey50",3),rep(NA,4)),
#          cex = 0.7,
#          pch = c(21,21,21,16,22,22,22),
#          pt.bg = c("#cc2a36", "#eb6841","#edc951","#ffffff00","#4d3434",
#                    "#4d343470","#4d343440"),
#          pt.cex = c(2,1,0.5,2,2,2,2))
#   if (!is.null(name)){
#     if(interactive()){
#       x <- locator()
#       if (!is.null(x)){
#         X <- sp::SpatialPoints(data.frame(x), proj4string = spdf@proj4string)
#         text(x = x$x, y = x$y, labels = sp::over(X, spdf)[,name], cex = 0.6,
#              adj = c(0,0))
#       }
#     } else {
#       text(x[x$col == "#cc2a36",c("X","Y")],labels = x[x$col == "#cc2a36", name],
#            cex = 0.6)
#     }
#   }
# }

