## ----, fig.show='hold', fig.width = 6, fig.height=6, results='hide'------
library(flow)

data(LoireAtlantique)

dom1<- flowDom(mat = MRE44,
               i = "DCRAN",
               j = "CODGEO",
               fij = "NBFLUX_C08_POP05P")

## ----, fig.show='hold', fig.width = 7, fig.height=7----------------------
head(dom1)

## ----, fig.show='hold', fig.width = 6, fig.height=6----------------------
plotflowDom(fdom = dom1, spdf = COM44, id = "INSEE_COM")
title("Residential Dominant Flows")

## ----, fig.show='hold', fig.width = 6, fig.height=6, results = 'hide'----
dom2 <- flowDom(mat = MRE44,
                i = "DCRAN",
                j = "CODGEO",
                fij = "NBFLUX_C08_POP05P", 
                k = 0.2)

## ----, fig.show='hold', fig.width = 6, fig.height=6----------------------
plotflowDom(fdom = dom2, spdf = COM44, id = "INSEE_COM")
title("Residential Dominant Flows (20% flows threshold)")

