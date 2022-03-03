#' @title Flow Selection from Origins
#' @name firstflows
#' @description Flow selection from origins.
#' @param mat A square matrix of flows.
#' @param method A method of flow selection, one of "nfirst", "xfirst" or "xsumfirst":
#' \itemize{
#' \item{nfirst selects the k first flows from origins,}
#' \item{xfirst selects flows greater than k,}
#' \item{xsumfirst selects as many flows as necessary for each origin so that their sum is at least equal to k.
#' If k is not reached for one origin, all its flows are selected.}
#' }
#' @param ties.method In case of equality with "nfirst" method, use "random" or "first" (see \link{rank}).
#' @param k Selection threshold.
#' @return A boolean matrix of selected flows.
#' @details As the output is a boolean matrix, use element-wise multiplication to get flows intensity.
#' @noRd
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
    matfinal[names(lfirst[i][]),lfirst[[i]][is.na(lfirst[[i]])==F]] <- 1
  }
  return(matfinal)
}


#' @title Flow Selection Based on Global Criteria
#' @name firstflowsg
#' @description Flow selection based on global criteria.
#' @param mat A square matrix of flows.
#' @param method A method of flow selection, one of "nfirst", "xfirst" or "xsumfirst":
#' \itemize{
#' \item{nfirst selects the k first flows of the matrix,}
#' \item{xfirst selects flows greater than k,}
#' \item{xsumfirst selects as many flows as necessary so that their sum is at least equal to k.}
#' }
#' @param ties.method In case of equality with "nfirst" method, use "random" or "first" (see \link{rank}).
#' @param k Selection threshold.
#' @return A boolean matrix of selected flows.
#' @details As the output is a boolean matrix, use element-wise multiplication to get flows intensity.
#' @noRd
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
#' @description Dominant flows selection.
#' @param mat A square matrix of flows.
#' @param w A vector of units weigths (sum of incoming flows, sum of outgoing flows...).
#' @param k A threshold (see 'Details').
#' @return A boolean matrix of selected flows.
#' @details This function selects which flow (fij or fji) must be kept.
#' If the ratio weight of destination (wj) / weight of origin (wi) is greater
#' than k, then fij is selected and fji is not.
#' This function can perform the second criterion of the Nystuen &
#' Dacey's dominants flows analysis.\cr
#' As the output is a boolean matrix, use element-wise multiplication to get flows intensity.
#' @references J. Nystuen & M. Dacey, 1961, "A Graph Theory Interpretation of Nodal Regions",
#' \emph{Papers and Proceedings of the Regional Science Association}, 7:29-42.
#' @examples
#' # Import data
#' data(nav)
#' myflows <- prepflows(mat = nav, i = "i", j = "j", fij = "fij")
#'
#' # Remove the matrix diagonal
#' diag(myflows) <- 0
#'
#' # Select the dominant flows (incoming flows criterion)
#' flowSel <- domflows(mat = myflows, w = colSums(myflows), k = 1)
#' statmat(mat = myflows * flowSel, output = "none")
#' @noRd
domflows <- function(mat, w, k){
  # list of i, j selected
  matfinal <- mat
  matfinal[] <- 0
  for (i in 1:dim(mat)[1]){
    for (j in 1:dim(mat)[2]){
      if (w[i] > 0){
        if ((w[j]/w[i]) > k){
          matfinal[i,j] <- 1
        }
      }
    }
  }
  matfinal[mat == 0] <- 0
  return(matfinal)
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


#' @title LegendPropLines
#' @name LegendPropLines
#' @import graphics
#' @noRd
LegendPropLines<- function(pos = "topleft", legTitle = "Title of the legend", legTitleCex = 0.8,
                           legValuesCex = 0.6, varvect, sizevect, col="red", frame=FALSE, round=0){


  positions <- c("bottomleft", "topleft", "topright", "bottomright",
                 "left", "right", "top", "bottom", "middle")
  if(pos %in% positions){

    # extent
    x1 <- par()$usr[1]
    x2 <- par()$usr[2]
    y1 <- par()$usr[3]
    y2 <- par()$usr[4]
    xextent <- x2 - x1
    yextent <- y2 - y1

    # variables internes
    paramsize1 <- 25
    paramsize2 <- 40
    width <- (x2 - x1) / paramsize1
    height <- width /1.5
    delta1 <- min((y2 - y1) / paramsize2, (x2 - x1) / paramsize2) # Gros eccart entre les objets
    delta2 <- (min((y2 - y1) / paramsize2, (x2 - x1) / paramsize2))/2 # Petit eccart entre les objets


    rValmax <- max(varvect,na.rm = TRUE)
    rValmin <- min(varvect,na.rm = TRUE)
    rValextent <- rValmax - rValmin
    rLegmax <- max(sizevect,na.rm = TRUE)
    rLegmin <- min(sizevect,na.rm = TRUE)
    rLegextent <- rLegmax - rLegmin

    rVal <- c(rValmax,rValmax - rValextent/3 , rValmax - 2*(rValextent/3),rValmin)
    rLeg <- c(rLegmax,rLegmax - rLegextent/3 , rLegmax - 2*(rLegextent/3),rLegmin)
    rVal <- round(rVal,round)

    # xsize & ysize

    longVal <- rVal[strwidth(rVal,cex=legValuesCex)==max(strwidth(rVal,cex=legValuesCex))][1]
    #if(!is.null(breakval)){if (strwidth(paste(">=",breakval),cex=legValuesCex)>strwidth(longVal,cex=legValuesCex)){longVal <- paste(">=",breakval)}}
    legend_xsize <- max(width+ strwidth(longVal,cex=legValuesCex)-delta2,strwidth(legTitle,cex = legTitleCex)-delta1)

    legend_ysize <-8*delta1 + strheight(legTitle,cex = legTitleCex)

    # Position
    if (pos == "bottomleft") {xref <- x1 + delta1 ; yref <- y1 + delta1}
    if (pos == "topleft") {xref <- x1 + delta1 ; yref <- y2 - 2*delta1 - legend_ysize}
    if (pos == "topright") {xref <- x2 - 2*delta1 - legend_xsize ; yref <- y2 -2*delta1 - legend_ysize}
    if (pos == "bottomright") {xref <- x2 - 2*delta1 - legend_xsize ; yref <- y1 + delta1}
    if (pos == "left") {xref <- x1 + delta1 ; yref <- (y1+y2)/2-legend_ysize/2 - delta2}
    if (pos == "right") {xref <- x2 - 2*delta1 - legend_xsize ; yref <- (y1+y2)/2-legend_ysize/2 - delta2}
    if (pos == "top") {xref <- (x1+x2)/2 - legend_xsize/2 ; yref <- y2 - 2*delta1 - legend_ysize}
    if (pos == "bottom") {xref <- (x1+x2)/2 - legend_xsize/2 ; yref <- y1 + delta1}
    if (pos == "middle") { xref <- (x1+x2)/2 - legend_xsize/2 ; yref <- (y1+y2)/2-legend_ysize/2 - delta2}


    # Frame
    if (frame==TRUE){
      rect(xref-delta1, yref-delta1, xref+legend_xsize + delta1*2, yref+legend_ysize + delta1 *2, border = "black",  col="white")
    }

    mycol <- col

    jump <- delta1
    for(i in 4:1){

      if (rLeg[i] < 0.2){rLeg[i] <- 0.2} # TAILLE DES LIGNE MINIMALES (A METTRE AUSSI SUR LES CARTES)

      segments(xref, yref + jump, xref + width, yref + jump, col=mycol, lwd=rLeg[i],lend=1)
      text(xref + width + delta2 ,y= yref + jump,rVal[i],adj=c(0,0.5),cex=legValuesCex)
      jump <- jump + 2*delta1 # ICI AMELIORER
    }
    text(x=xref ,y=yref + 9*delta1,legTitle,adj=c(0,0),cex=legTitleCex)
  }
}
