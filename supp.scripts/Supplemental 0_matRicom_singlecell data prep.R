### creation of single cell data objects ###

# PBMC data are from https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
dr <- getwd()
pbmc <- make.object("pbmc",dr)

# HNSCC data are from NicheNet
# these data are not in 10X format, hence the make.object function cannot be used. See below.

# CRC data are from https://www.10xgenomics.com/resources/datasets
dr <- getwd() #make sure the directory and sub-directory are structured as explained in the function help page
CRC <- make.object.spatial(dr,"filtered_feature_bc_matrix.h5")

## ALTERNATIVELY, STEP BY STEP ##

#pbmc object from Seurat vignette
library(dplyr)
library(Seurat)
library(patchwork)
pbmc.data <- Read10X(data.dir = "../data/pbmc3k/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc,"pbmc.RDS")

#hnscc object from Nichenetr vignette
hnscc <- readRDS(url("https://zenodo.org/record/3260758/files/hnscc_expression.rds"))
expression <- hnscc$expression
sample_info <- hnscc$sample_info 
tumors_remove <- c("HN10","HN","HN12", "HN13", "HN24", "HN7", "HN8","HN23")
CAF_ids <- sample_info %>% filter(`Lymph node` == 0) %>% filter((tumor %in% tumors_remove == FALSE)) %>% filter(`non-cancer cell type` == "CAF") %>% .$cell
malignant_ids <- sample_info %>% filter(`Lymph node` == 0) %>% filter(`classified  as cancer cell` == 1) %>% filter((tumor %in% tumors_remove == FALSE)) %>% .$cell
expressed_genes_CAFs <- expression[CAF_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
expressed_genes_malignant <- expression[malignant_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
#cutting to HN5 tumor for brevity, as for the original Nichenetr vignette
sample_info <- subset(sample_info,sample_info$tumor=="HN5")
expression <- expression[rownames(expression)%in%sample_info$cell,]
#create a Seurat object
hnscc <- CreateSeuratObject(counts = t(expression), meta.data = sample_info)
hnscc[["percent.mt"]] <- PercentageFeatureSet(hnscc, pattern = "^MT-")
hnscc <- FindVariableFeatures(hnscc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(hnscc)
hnscc <- ScaleData(hnscc, features = all.genes)
hnscc <- RunPCA(hnscc, features = VariableFeatures(object = hnscc))
hnscc <- FindNeighbors(hnscc, dims = 1:10)
hnscc <- FindClusters(hnscc, resolution = 0.5)
hnscc <- RunUMAP(hnscc, dims = 1:10)
hnscc@meta.data$classified..as.cancer.cell <- sample_info$`classified  as cancer cell`
hnscc@meta.data$classified.as.non.cancer.cells <- sample_info$`classified as non-cancer cells`
hnscc@meta.data$non.cancer.cell.type <- sample_info$`non-cancer cell type`
hnscc@meta.data$non.cancer.cell.type <- ifelse(hnscc@meta.data$non.cancer.cell.type==0,"cancer.cell",hnscc@meta.data$non.cancer.cell.type)
saveRDS(hnscc,"hnscc.RDS")

#spatial CRC object from 10X genomics - prepared following the standard spatial Seurat workflow
crc <- Load10X_Spatial(data.dir = getwd())
crc <- SCTransform(crc, assay = "Spatial", verbose = FALSE)
crc <- RunPCA(crc, assay = "SCT", verbose = FALSE)
crc <- FindNeighbors(crc, reduction = "pca", dims = 1:30)
crc <- FindClusters(crc, verbose = FALSE)
crc <- RunUMAP(crc, reduction = "pca", dims = 1:30)
saveRDS(crc,"CRC_visio10X_slide1.RDS")

