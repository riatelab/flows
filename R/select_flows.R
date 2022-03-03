#' @title Flow selection
#' @name select_flows
#' @description Flow selection from origins.
#' @param mat A square matrix of flows.
#' @param method A method of flow selection, one of "dominant", "nfirst",
#' "xfirst" or "xsumfirst":
#' \itemize{
#' \item{dominant selects the dominant flows (see Details)}
#' \item{nfirst selects the k first flows from origins,}
#' \item{xfirst selects flows greater than k,}
#' \item{xsumfirst selects as many flows as necessary for each origin so that their sum is at least equal to k.
#' If k is not reached for one origin, all its flows are selected.}
#' }
#' @param ties In case of equality with "nfirst" method, use "random" or "first" (see \link{rank}).
#' @param global If TRUE flows selections is done at the matrix scale.
#' @param k Selection threshold for nfirst, xfirst and xsumfirst methods,
#' ratio for dominant method.
#' @param w A vector of units weigths (sum of incoming flows, sum of outgoing flows...).
#' @export
#' @return
#' A boolean matrix of selected flows.
#' Use element-wise multiplication to get flows intensity.
#' @details
#' If method = "dominant", select which flow (fij or fji) must be kept.
#' If the ratio weight of destination (wj) / weight of origin (wi) is greater
#' than k, then fij is selected and fji is not.
#' This function can perform the second criterion of the Nystuen &
#' Dacey's dominants flows analysis.
#' @examples
#' # Import data
#' nav <- read.csv(system.file("csv/nav.csv", package = "flows"))
#' # Prepare data
#' mat <- prepare_mat(x = nav, i = "i", j = "j", fij = "fij")
#' # remove diagonal
#' diag(mat) <- 0
#'
#' # Select the first flow from each origin
#' res <- select_flows(mat = mat, method = "nfirst", global = FALSE, k = 1)
#' rowSums(res)
#'
#' # Select the 5 first flows of the matrix (from origins)
#' res <- select_flows(mat = mat, method = "nfirst", global = TRUE, k = 5)
#' sum(rowSums(res))
#'
#' # Select the flows greater than 5000
#' res <- select_flows(mat = mat, method = "xfirst", k = 5000)
#' r <- mat * res
#' r[r>0]
#'
#' # Select as many flows as necessary for each origin so that their sum is at least equal to 500.
#' res <- select_flows(mat = mat, method = "xsumfirst", global = FALSE, k = 500)
#' r <- mat * res
#' rowSums(r)
#'
#' # Select as many flows in the matrix so that their sum is at least equal to 50000.
#' res <- select_flows(mat = mat, method = "xsumfirst", global = TRUE, k = 50000)
#' r <- mat * res
#' sum(rowSums(r))
#'
#' # Select dominant flows
#' m <- mat[1:5,1:5]
#' ws <- colSums(m)
#' res <- select_flows(mat = m, method = "dominant", k = 1, w = ws)
#' # 2nd element has a lower weight than 3rd element (ratio > 1)
#' ws[3] / ws[2]
#' # The flow from 2nd element to 3rd element is kept
#' res[2, 3]
#' # The flow from 3rd element to 2nd element is removed
#' res[3, 2]
select_flows <- function(mat, method = "nfirst", ties = "first", global = FALSE, k, w){
  if(method %in% c('nfirst', "xfirst", "xsumfirst")){
    if(global == TRUE){
      x <- firstflowsg(mat = mat, method = method,  k = k, ties.method = ties)
    }
    if(global == FALSE){
      x <- firstflows(mat = mat, method = method,  k = k, ties.method = ties)
    }
  }
  if (method == "dominant"){
    x <- domflows(mat = mat, w = w, k = k)
  }
  return(x)
}
