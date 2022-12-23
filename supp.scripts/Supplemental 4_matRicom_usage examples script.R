### Usage example script ###
# for info on the datasets, see Supplementary script 0 #
# All data are available online at https://doi.org/10.5281/zenodo.7385244 #

library(matRicom)

#example 1
pbmc <- readRDS("pbmc.RDS")
results <- matricom(pbmc, "new.ids")

flowgraph(results) #equal to flowgraph(results,color.by = "sender",targets = T)
flowgraph(results, color.by = "ecm")
flowgraph(results, color.by = "receiver")

annograph(results) #equal to annograph(results,select = "sender")
annograph(results, select = "receiver")

#example 2
results <- matricom(pbmc, "new.ids", explainable = F)
flowgraph(results, simplify = T, targets = F)
annograph(results)

#example 3
hnscc <- readRDS("hnscc.RDS")
results <- matricom(hnscc, "non.cancer.cell.type")
flowgraph(results)
flowgraph(results, simplify = T, targets = F, color.by = "ecm")
annograph(results)

#example 4
results <- matricom(hnscc, "non.cancer.cell.type", target = "CAF", explainable = T)
flowgraph(results, color.by = "sender")
flowgraph(results, color.by = "ecm", simplify = F)
annograph(results, select = "receiver")

#example 5
#note that the compgraph function is highly experimental and results aren't guaranteed to make sense
h1 <- matricom(hnscc, "non.cancer.cell.type", target = "CAF", explainable = F)
h2 <- matricom(hnscc, "non.cancer.cell.type", target = "Macrophage", explainable = F)
h3 <- matricom(hnscc, "non.cancer.cell.type", target = "Endothelial", explainable = F)
l <- list(caf=h1,
          mo=h2)
compgraph(l, targets = F)
compgraph(l, targets = T, only.unique = F)
l <- list(caf=h1,
          endothelial=h3)
compgraph(l, only.unique = F)
compgraph(l, targets = T, only.unique = F)

#example 6
results <- matricom(hnscc, "non.cancer.cell.type", target = "CAF", explainable = T)
netgraph(results) #random picks for ECM and receiving population
netgraph(results, receiving.population = "CAF", ECM.component = "FN1")
netgraph(results, label.size = 0.5)
netgraph(results, label.size = 0.5, simplified = F)
netgraph(results, interactive = TRUE)

#example 7
CRC <- readRDS("CRC_visio10X_slide1.RDS") 
results <- matricom.spatial(CRC,"seurat_clusters")
spatplot(results,CRC,"seurat_clusters")
spatplot(results,CRC,"seurat_clusters","COL1A1_CD44_MT-CO1_ZBTB20")
spatplot(results,CRC,"seurat_clusters","COL1A1_CD44_MT-CO1_ZBTB20",singles=T)
