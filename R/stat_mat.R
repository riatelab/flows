#' @title Descriptive statistics on flow matrix
#' @name stat_mat
#' @description This function provides various indicators and graphical outputs
#' on a flow matrix.
#' @param mat A square matrix of flows.
#' @param output Graphical output. Choices are "all" for all graphics,
#' "none" to avoid any graphical output, "degree" for degree distribution, "wdegree" for
#' weighted degree distribution, "lorenz" for Lorenz curve of link weights and
#' "boxplot" for boxplot of link weights (see 'Details').
#' @param verbose A boolean, if TRUE, returns statistics in the console.
#' @return  The function returns a list of statistics and may plot graphics.
#' \itemize{
#' \item{nblinks: number of cells with values > 0}
#' \item{density: number of links divided by number of possible links (also called gamma index by geographers), loops excluded}
#' \item{connectcomp: number of connected components (isolates included,
#' weakly connected: use of \code{\link{clusters}} where mode = "weak")}
#' \item{connectcompx: number of connected components (isolates deleted,
#' weakly connected: use of \code{\link{clusters}} where mode = "weak")}
#' \item{sizecomp: a data.frame of connected components: size
#' and sum of flows per component (isolates included)}
#' \item{compocomp: a data.frame of connected components giving membership of units (isolates included)}
#' \item{degrees: a data.frame of nodes degrees and weighted degrees}
#' \item{sumflows: sum of flows}
#' \item{min: minimum flow }
#' \item{Q1: first quartile of flows}
#' \item{median: median flow}
#' \item{Q3: third quartile of flows}
#' \item{max: maximum flow}
#' \item{mean: mean flow}
#' \item{sd: standart deviation of flows}}
#' @details Graphical ouputs concern outdegrees by default. If the matrix is
#' transposed, outputs concern indegrees.
#' @import graphics
#' @import stats
#' @examples
#' # Import data
#' nav <- read.csv(system.file("csv/nav.csv", package = "flows"))
#' myflows <- prepare_mat(x = nav, i = "i", j = "j", fij = "fij")
#'
#' # Get statistics and graphs about the matrix
#' mystats <- stat_mat(mat = myflows, output = "all", verbose = TRUE)
#'
#' # Size of connected components
#' mystats$sizecomp
#'
#' # Sum of flows
#' mystats$sumflows
#'
#' # Plot Lorenz curve only
#' stat_mat(mat = myflows, output = "lorenz", verbose = FALSE)
#'
#' # Statistics only
#' mystats <- stat_mat(mat = myflows, output = "none", verbose = FALSE)
#' str(mystats)
#' @export
stat_mat <- function(mat, output = "all", verbose = TRUE){
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

  #prep graph
  deg <- rowSums(matbool)
  deg2 <- rowSums(mat)
  # df output
  degdf <- data.frame(id=names(deg), degree = deg, wdegree = deg2, stringsAsFactors = FALSE)
  deg <- deg[deg>0]
  deg2 <- deg2[deg2>0]

  if(output=="degree"){
    plot(deg[order(deg, decreasing = TRUE)], type = "l", log = "xy",
         xlab = "rank (log)", ylab = "size (log nb. flows)")
    title("rank - size")
  }
  if(output=="wdegree"){
    plot(deg2[order(deg2, decreasing = TRUE)], type = "l", log = "xy",
         xlab = "rank (log)", ylab = "size (log flow intensity)")
    title("rank - size (weighted)")
  }
  if(output=="lorenz"){
    plot( y = vmatcs, x = seq(0,100,length.out = length(vmatcs)), type = "l",
          xlim = c(0,100), ylim = c(0,100),
          xlab = "cum. nb. flows", ylab = "cum. intensity of flows")
    title ("Lorenz Curve")
  }
  if(output=="boxplot"){
    boxplot(as.vector(mat[mat>0]), log = "y")
    title("Boxplot")
  }
  if(output=="all"){
    ## graphic outputs
    old.par <- par (mfrow = c(2,2))
    ## rank-size link
    plot(deg[order(deg, decreasing = TRUE)], type = "l", log = "xy",
         xlab = "rank (log)", ylab = "size (log nb. flows)")
    title("rank - size")
    ## rank size flow
    plot(deg2[order(deg2, decreasing = TRUE)], type = "l", log = "xy",
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
  }

  ## Connected components of a graph
  g <- igraph::graph.adjacency(adjmatrix = mat, mode = "directed", weighted = TRUE)
  clustg <- igraph::clusters(graph = g, mode = "weak")
  connectcomp <- clustg$no
  connectcompx <- length(clustg$csize[clustg$csize>1])
  compocomp <-  data.frame(id = V(g)$name, idcomp = clustg$membership, stringsAsFactors = FALSE)
  compocompw <- merge(compocomp, degdf, by = "id")
  compw <- aggregate(x = compocompw$wdegree,by = list(compocompw$idcomp),
                     FUN = sum)
  sizecomp <- data.frame(idcomp = seq(1, length(clustg$csize)),
                         sizecomp = clustg$csize, wcomp = compw$x)

  if (verbose == TRUE){
    ## stat cat
    cat('matrix dimension:', matdim, "X", matdim,"\n" )
    cat('nb. links:', nbcellfull, "\n" )
    cat('density:', nbcellfull/(nbcell - matdim), "\n" )
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
  }
  ## stat list
  matstat <- list(matdim = dim(mat),
                  nblinks = nbcellfull,
                  density = nbcellfull/(nbcell - matdim),
                  connectcomp = connectcomp,
                  connectcompx = connectcompx,
                  sizecomp = sizecomp,
                  compocomp = compocomp,
                  degrees = degdf,
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
