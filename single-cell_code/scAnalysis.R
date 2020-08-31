#you may need to install these packages:
#install.packages("dplyr")
#install.packages("Seurat")
#install.packages("patchwork")

#you also need to download the data from 
#https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz,
#extract it somewhere and change directory to that folder
#e.g.: 
setwd("C:/Work/Teaching/KMG-060/KMG060-Systems-Biology-course/single-cell_code/data")

library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data). 
# - We already filter out the empty droplets here, by setting mininal features, and remove genes only expressed in a couple of cells or less
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
#calculate mitochondrial content
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#show the quality metrics against each other
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#we see that a few cells has a very high mitochondrial content

#filter bad cells
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#normalize and log transform
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

#find highly variable genes

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels (selects 2000 genes)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

#scaling
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

#show which genes that are primarily included in PC1 and PC2:
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

#display the 2 first PCs in 2D:
DimPlot(pbmc, reduction = "pca")

#get an overview of how the cells are separated by the genes in the PCs
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

#Try to figure out how many of the PCs to use - the ones that explain little variance will mostly add noise
#We will not run Jackstraw, since it takes 5 minutes to run
pbmc <- JackStraw(pbmc, num.replicate = 100) #this takes a while to run
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

JackStrawPlot(pbmc, dims = 1:15)

#We instead use something faster and simpler, called an elbow plot
ElbowPlot(pbmc)


#cluster:
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)


#UMAP for visaulization
pbmc <- RunUMAP(pbmc, dims = 1:10)

#visualize clusters - the umap may change between runs, so don't worry about that
DimPlot(pbmc, reduction = "umap")

# find all markers of cluster 1 to be able to figure out the cell type
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

#we skip some additional searches for markers here

#Now, plot which cells that express some of the markers
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                               "CD8A"))

new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
