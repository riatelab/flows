#' @title Comparison of two matrices
#' @name compare_mat
#' @description Compares two matrices of same dimension, with same column and
#' row names.
#' @param mat1 A square matrix of flows.
#' @param mat2 A square matrix of flows.
#' @param digits	An integer indicating the number of decimal places to be used
#' when printing the data.frame in the console (see \link{round}).
#' @return A data.frame that provides statistics on differences
#' between mat1 and mat2: absdiff are the
#' absolute differences and reldiff are the relative differences (in percent).
#' @examples
#' # # Import data
#' nav <- read.csv(system.file("csv/nav.csv", package = "flows"))
#' mat <- prepare_mat(x = nav, i = "i", j = "j", fij = "fij")
#' # Remove the matrix diagonal
#' diag(mat) <- 0
#'
#' # Select the first flows
#' flowSel1 <- select_flows(mat = mat, method = "nfirst", k = 1)
#'
#' # Select flows greater than 2000
#' flowSel2 <- select_flows(mat = mat, method = "xfirst", k = 2000)
#'
#' # Combine selections
#' flowSel <- mat * flowSel1 * flowSel2
#'
#' # Compare flow matrices
#' compare_mat(mat1 = mat, mat2 = flowSel, digits = 1)
#' @export
compare_mat <- function(mat1, mat2, digits = 0){
  x1 <- stat_mat(mat1, output = "none", verbose = FALSE)
  x2 <- stat_mat(mat2, output = "none", verbose = FALSE)
  compdf <- data.frame(mat1= c(x1$nblinks, x1$sumflows, x1$connectcompx,
                               x1$min, x1$Q1, x1$median, x1$Q3,
                               x1$max, x1$mean, x1$sd),
                       mat2= c(x2$nblinks, x2$sumflows,x2$connectcompx,
                               x2$min, x2$Q1, x2$median, x2$Q3,
                               x2$max, x2$mean, x2$sd),
                       row.names = c("nblinks","sumflows", "connectcompx",
                                     "min", "Q1", "median", "Q3",
                                     "max", "mean", "sd"))

  compdf$absdiff <- abs(compdf$mat1-compdf$mat2)
  compdf[4:10,"absdiff"] <- NA
  compdf$reldiff <- abs(compdf$mat1-compdf$mat2) / compdf$mat1 * 100
  compdf[3:10,"reldiff"] <- NA

  print(round(compdf, digits))
  return(invisible(compdf))
}
