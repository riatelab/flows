## ----, fig.show='hold', fig.width = 7, fig.height=7----------------------
library(flow)

data(lot)

dom1<- flowDom(mat = MRE46,
               i = "DCRAN",
               j = "CODGEO",
               fij = "NBFLUX_C08_POP05P")

head(dom1)

plotflowDom(fdom = dom1, spdf = COM46, id = "INSEE_COM")
title("Residential Dominant Flows")

## ----, fig.show='hold', fig.width = 7, fig.height=7----------------------

dom2 <- flowDom(mat = MRE46,
                i = "DCRAN",
                j = "CODGEO",
                fij = "NBFLUX_C08_POP05P",
                k = 0.5)

plotflowDom(fdom = dom2, spdf = COM46, id = "INSEE_COM")
title("Residential Dominant Flows (50% flows threshold)")

## ----, fig.show='hold', fig.width = 7, fig.height=7----------------------

plotflowDom(fdom = dom2, spdf = COM46, id = "INSEE_COM", name = "NOM_COM")
title("Residential Dominant Flows (50% flows threshold)")

