% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stat_mat.R
\name{stat_mat}
\alias{stat_mat}
\title{Descriptive statistics on flow matrix}
\usage{
stat_mat(mat, output = "all", verbose = TRUE)
}
\arguments{
\item{mat}{A square matrix of flows.}

\item{output}{Graphical output. Choices are "all" for all graphics,
"none" to avoid any graphical output, "degree" for degree distribution, "wdegree" for
weighted degree distribution, "lorenz" for Lorenz curve of link weights and
"boxplot" for boxplot of link weights (see 'Details').}

\item{verbose}{A boolean, if TRUE, returns statistics in the console.}
}
\value{
The function returns a list of statistics and may plot graphics.
\itemize{
\item{nblinks: number of cells with values > 0}
\item{density: number of links divided by number of possible links (also called gamma index by geographers), loops excluded}
\item{connectcomp: number of connected components (isolates included,
weakly connected: use of \code{\link{clusters}} where mode = "weak")}
\item{connectcompx: number of connected components (isolates deleted,
weakly connected: use of \code{\link{clusters}} where mode = "weak")}
\item{sizecomp: a data.frame of connected components: size
and sum of flows per component (isolates included)}
\item{compocomp: a data.frame of connected components giving membership of units (isolates included)}
\item{degrees: a data.frame of nodes degrees and weighted degrees}
\item{sumflows: sum of flows}
\item{min: minimum flow }
\item{Q1: first quartile of flows}
\item{median: median flow}
\item{Q3: third quartile of flows}
\item{max: maximum flow}
\item{mean: mean flow}
\item{sd: standart deviation of flows}}
}
\description{
This function provides various indicators and graphical outputs
on a flow matrix.
}
\details{
Graphical ouputs concern outdegrees by default. If the matrix is
transposed, outputs concern indegrees.
}
\examples{
# Import data
nav <- read.csv(system.file("csv/nav.csv", package = "flows"))
myflows <- prepare_mat(x = nav, i = "i", j = "j", fij = "fij")

# Get statistics and graphs about the matrix
mystats <- stat_mat(mat = myflows, output = "all", verbose = TRUE)

# Size of connected components
mystats$sizecomp

# Sum of flows
mystats$sumflows

# Plot Lorenz curve only
stat_mat(mat = myflows, output = "lorenz", verbose = FALSE)

# Statistics only
mystats <- stat_mat(mat = myflows, output = "none", verbose = FALSE)
str(mystats)
}
