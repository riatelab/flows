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
firstflows <- function(mat, method = "nfirst", ties.method = "first", k) {
  # list of i, j selected
  lfirst <- apply(mat, 1, get(method), k = k, ties.method = ties.method)
  # if only one selected
  if (is.null(dim(lfirst))) {
    lfirst <- as.list((lfirst))
  }
  # control class output
  if (is.matrix(lfirst)) {
    lfirst <- as.list(as.data.frame(lfirst, stringsAsFactors = FALSE))
  }
  matfinal <- mat
  matfinal[] <- 0
  # control 0 selection
  if (length(lfirst) < 1) {
    return(matfinal)
  }
  for (i in 1:nrow(matfinal)) {
    matfinal[names(lfirst[i][]), lfirst[[i]][is.na(lfirst[[i]]) == FALSE]] <- 1
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
firstflowsg <- function(mat, method = "nfirst", k, ties.method = "first") {
  matfinal <- mat
  matfinal[] <- 0
  if (method == "nfirst") {
    matfinal[rank(mat, ties.method = ties.method) > ((dim(mat)[1] * dim(mat)[2]) - k)] <- 1
  }
  if (method == "xfirst") {
    matfinal[mat >= k] <- 1
  }
  if (method == "xsumfirst") {
    matv <- as.vector(mat)
    names(matv) <- 1:length(matv)
    matvo <- matv[order(matv, decreasing = TRUE)]
    matvo <- cumsum(matvo)
    nbgood <- (length(matvo[matvo < k]) + 1)
    matvo[] <- c(rep(1, nbgood), rep(0, (length(matvo) - nbgood)))
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
domflows <- function(mat, w, k) {
  # list of i, j selected
  matfinal <- mat
  matfinal[] <- 0
  for (i in 1:dim(mat)[1]) {
    for (j in 1:dim(mat)[2]) {
      if (w[i] > 0) {
        if ((w[j] / w[i]) > k) {
          matfinal[i, j] <- 1
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
nfirst <- function(x, k, ties.method) {
  x <- x[x > 0]
  if (length(x) > 0) {
    if (length(x) > k) {
      x <- x[rank(x, ties.method = ties.method) > (length(x) - k)]
      x <- names(x)
    } else {
      x <- names(x)
    }
  }
  return(x)
}

#' @title xfirst
#' @name xfirst
#' @noRd
xfirst <- function(x, k, ties.method) {
  x <- x[x > 0]
  if (length(x) > 0) {
    x <- x[x > k]
    x <- names(x)
  } else {
    x <- names(x)
  }
  return(x)
}

#' @title xsumfirst
#' @name xsumfirst
#' @noRd
xsumfirst <- function(x, k, ties.method) {
  x <- x[x > 0]
  x <- x[order(x, decreasing = TRUE)]
  x <- cumsum(x = x)
  if (length(x) > 0) {
    if (x[1] >= k) {
      x <- names(x[1])
    } else {
      # at least k
      x <- x[1:(length(x[(x <= k) == TRUE]) + 1)]
      # less than k
      # x <- x[(x <= k)]
      x <- names(x)
    }
  } else {
    x <- names(x)
  }
  # pb return error if k > sum(mat) or k > sum(mat/w)
  return(x)
}

