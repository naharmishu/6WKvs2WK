## 6 weeks HDM data vs 2 weeks HDM data and dog allergen data
## Single-cell RNA_Seq datasets in R using Seurat
# JEBUNNAHAR MISHU

## Loading the libraries
library(BiocManager)
library(devtools)
library(Seurat)
library(SeuratObject)
library(remotes)
library(targets)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(Matrix)
library(pheatmap)
library(harmony)


# Get or read the file
base_path <- "/Users/cbm737/Documents/development"
# Set the directory
setwd(base_path)

Dirs <- list.dirs(path = base_path, full.names = FALSE, recursive = FALSE)
Dirs
list.files(Dirs)


for (x in Dirs) {
  name <- gsub('_filtered_feature_bc_matrix','', x)
  
  # Construct file paths
  mtx_path <- file.path(base_path, x, "matrix.mtx.gz")
  features_path <- file.path(base_path, x, "features.tsv.gz")
  barcodes_path <- file.path(base_path, x, "barcodes.tsv.gz")
  
  # Read the data
  DAHDM_count <- ReadMtx(mtx = mtx_path, 
                         features = features_path, 
                         cells = barcodes_path)
  # create seurat object
  assign(name, CreateSeuratObject(count = DAHDM_count))
  
}

ls()

### merge datasets
HDMDA_Obj <- merge(HDM6WK, y = c(HDM2WK, DA2WK), add.cell.ids = c("HDM6_6WK", "HDM2_2WK", "DA2_2WK"), project = "HDMvsDA") 
str(HDMDA_Obj)
View(HDMDA_Obj@meta.data)
head(HDMDA_Obj@meta.data)
HDMDA_Obj

# create the column
HDMDA_Obj$Sample <- rownames(HDMDA_Obj@meta.data)

# splite the column
HDMDA_Obj@meta.data <- separate(HDMDA_Obj@meta.data, col = 'Sample', into = c('Id', 'Time', 'Barcode'), sep = '_')
unique(HDMDA_Obj@meta.data$Id)

# Joinlayers
HDMDA_Obj <- JoinLayers(HDMDA_Obj)

# percentage mitocondrial 
HDMDA_Obj[["PercentMT"]] <- PercentageFeatureSet(HDMDA_Obj, pattern = "^mt-")

## Removing Tcr  
HDMDA_Obj <- HDMDA_Obj[!grepl('^Tr[abdg][vjc]', rownames(HDMDA_Obj))] # mouse
HDMDA_Obj

## Visualize data before filtering as a violin plot
VlnPlot(HDMDA_Obj, features = c("nCount_RNA", "nFeature_RNA","PercentMT"), ncol = 3)

#filtering
HDMDA_Obj_filtered <- subset(HDMDA_Obj, subset = nCount_RNA < 20000 & nFeature_RNA > 200 & PercentMT < 15)
HDMDA_Obj_filtered

## Visualize data before normalization as a violin plot
VlnPlot(HDMDA_Obj_filtered, features = c("nCount_RNA", "nFeature_RNA","PercentMT"), ncol = 3)
VlnPlot(HDMDA_Obj_filtered, features = c("nCount_RNA", "nFeature_RNA","PercentMT"), ncol = 3, pt.size = 0)

# Feature Scatter is used to visualize feature to feature relations
FeatureScatter(HDMDA_Obj_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = lm)

# Normalize the data
## Find the variables features
## Scale the data to avoid technical noise (batch effect) or biological sources ( different cell cycle)
## perform linear dimension reduction
# Find neighbors, clusters and Run UMAP to visualization
HDMDA_Obj_Norm <- NormalizeData(object = HDMDA_Obj_filtered) %>%
                   FindVariableFeatures() %>% 
                    ScaleData() %>%
                     RunPCA() %>%
                      FindNeighbors(dims = 1:17) %>%
                       FindClusters(resolution = c(0.3, 0.5, 0.8)) %>%
                        RunUMAP(dims = 1:17)

str(HDMDA_Obj_Norm)
view(HDMDA_Obj_Norm@meta.data)


top_Vgenes <- head(VariableFeatures(HDMDA_Obj_Norm), 20)
plot1 <- VariableFeaturePlot(HDMDA_Obj_Norm)
plot2 <- LabelPoints(plot = plot1, points = top_Vgenes, repel = TRUE)
plot1 + plot2

ElbowPlot(HDMDA_Obj_Norm)

PCHeatmap(HDMDA_Obj_Norm, dim = 1:6, cells = 500, balanced =TRUE, ncol = 3)

UMAPPlot(HDMDA_Obj_Norm)

DimPlot(HDMDA_Obj_Norm, reduction = "umap", group.by = "RNA_snn_res.0.5", label = TRUE, label.box = TRUE)+ labs(title = "Without Integration")
DimPlot(HDMDA_Obj_Norm, reduction = "umap", group.by = "RNA_snn_res.0.5", split.by = "Id", label = TRUE, label.box = TRUE)+ labs(title = "Without Integration")

### INTEGRATION WITH HARMONY

HDMDA_Obj_Harm <- NormalizeData(object = HDMDA_Obj_filtered) %>%
                   FindVariableFeatures() %>% 
                    ScaleData() %>%
                     RunPCA() %>%
                      RunHarmony(group.by.vars = "Id") %>%
                       FindNeighbors(reduction = "harmony", dims = 1:17) %>%
                        FindClusters(resolution = c(0.3, 0.5, 0.8)) %>%
                         RunUMAP(reduction = "harmony", dims = 1:17)

str(HDMDA_Obj_Harm)
view(HDMDA_Obj_Harm@meta.data)


DimPlot(HDMDA_Obj_Harm, reduction = "umap", group.by = "RNA_snn_res.0.5", label = TRUE, label.box = TRUE)+ labs(title = "With Integration")
DimPlot(HDMDA_Obj_Harm, reduction = "umap", group.by = "RNA_snn_res.0.5", split.by = "Id", label = TRUE, label.box = TRUE)+ labs(title = "With Integration")

# Annotate cell cluster
Cl_marker <- FindAllMarkers(HDMDA_Obj_Harm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(Cl_marker)
dim(Cl_marker)


#View top markers for the specific cluster
Top10 <- Cl_marker %>% group_by(cluster) %>% top_n(n=10, wt = avg_log2FC)
view(Top10)

DoHeatmap(HDMDA_Obj_Harm, features = Top10$gene) + NoLegend()

write.csv(Top10, file = "Top10.csv")
