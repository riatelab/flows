## -----------------------------------------------------------------------------
library(flows)
# Import data
nav <- read.csv(system.file("csv/nav.csv", package = "flows"))
head(nav, 4)
# Prepare data
mat <- prepare_mat(x = nav, i = "i", j = "j", fij = "fij")
mat[1:4,1:4]

## ---- fig.width=5, fig.height=5-----------------------------------------------
# Import data
nav <- read.csv(system.file("csv/nav.csv", package = "flows"))
# Prepare data
mat <- prepare_mat(x = nav, i = "i", j = "j", fij = "fij")

# Get statistics about the matrix
stat_mat(mat = mat, output = "none", verbose = TRUE)

# Plot Lorenz curve only
stat_mat(mat = mat, output = "lorenz", verbose = FALSE)

## ---- fig.width=7, fig.height=7-----------------------------------------------
# Graphics only
stat_mat(mat = mat, output = "all", verbose = FALSE)

# Statistics only
mystats <- stat_mat(mat = mat, output = "none", verbose = FALSE)
str(mystats)
# Sum of flows
mystats$sumflows


## ---- fig.height = 4, fig.width=4---------------------------------------------
# Import data
nav <- read.csv(system.file("csv/nav.csv", package = "flows"))
# Prepare data
mat <- prepare_mat(x = nav, i = "i", j = "j", fij = "fij")

# Remove the matrix diagonal
diag(mat) <- 0

# Selection of flows > 500
mat_sel_1 <- select_flows(mat = mat, method = "xfirst", k = 500)

# Selection of flows > 1000
mat_sel_2 <- select_flows(mat = mat, method = "xfirst", k = 1000)

# Compare initial matrix and selected matrices
compare_mat(mat1 = mat, mat2 = mat * mat_sel_1, digits = 0)
compare_mat(mat1 = mat, mat2 = mat * mat_sel_2, digits = 0)

## ---- fig.height = 6, fig.width=6---------------------------------------------
# Import data
nav <- read.csv(system.file("csv/nav.csv", package = "flows"))
# Prepare data
mat <- prepare_mat(x = nav, i = "i", j = "j", fij = "fij")

# Remove the matrix diagonal
diag(mat) <- 0

# Percentage of each outgoing flows
mat_p <- mat  / rowSums(mat) * 100

# Select flows that represent at least 20% of the sum of outgoing flows for 
# each urban area.
mat_p_sel <- select_flows(mat = mat_p, method = "xfirst", k = 20)

# Compare initial and selected matrices
compare_mat(mat1 = mat, mat2 = mat * mat_p_sel)


## -----------------------------------------------------------------------------
nav <- read.csv(system.file("csv/nav.csv", package = "flows"))
mat <- prepare_mat(x = nav, i = "i", j = "j", fij = "fij")
res <- nodal_flows(mat)

compare_mat(mat1 = mat, mat2 = res)

## ---- fig.height=7, fig.width=7, eval = TRUE----------------------------------
library(sf)
library(mapsf)
nav <- read.csv(system.file("csv/nav.csv", package = "flows"))
mat <- prepare_mat(x = nav, i = "i", j = "j", fij = "fij")
UA <- st_read(system.file("gpkg/GE.gpkg", package = "flows"), 
              layer = "urban_area", quiet = TRUE)
GE <- st_read(system.file("gpkg/GE.gpkg", package = "flows"), 
              layer = "region", quiet = TRUE)
mf_map(GE, col = "#cceae7", border = NA)
out <- map_nodal_flows(mat = mat, x = UA,
                       col_node = c('red', 'orange', 'yellow'),
                       col_flow = "grey30",
                       breaks = c(4,100,1000,2500,8655),
                       lwd = c(1,4,8,16), add = TRUE)
mf_title("Dominant Flows of Commuters")
mf_credits("INSEE, 2011")
mf_scale()


head(out$nodes[order(out$nodes$w, decreasing = TRUE), 2:3, drop = TRUE])


