getwd()
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SingleR)
library(celldex)
library(RColorBrewer)
#Loading the 10x Genomic data and creating a Seurat object
#In NCBI GEO data set (GSE213625)
#Download barcode.tsv, genes.tsv, matrix.mtx files
#Set seed is used a key to get the same output/ reproducible
set.seed(123)

#Need to set the data directory for control
data.dir <- "data path address"
raw_data <- Read10X(data.dir=data.dir)
list.files(data.dir)

#Create Seurat object for control 
seurat_obj1 <-CreateSeuratObject(
  counts = raw_data,
  project = 'Arabidopsis_c',
  min.cells = 3, # genes detected in >3 cells
  min.features = 300 #cell with >200 genes
)
seurat_obj1
head(colnames(seurat_obj1))
#Add metadata
seurat_obj1$sample <-'control'
seurat_obj1$condition <- 'control'
head(seurat_obj1@meta.data)

#Add mitochondrial Percentage data
seurat_obj1[["percent.mt"]] <- PercentageFeatureSet(
  seurat_obj1,
  pattern = "^ATMG"
)
#Add chloroplast Percentage data
seurat_obj1[["percent.chloro"]] <- PercentageFeatureSet(
  seurat_obj1,
  pattern = "^ATCG"
)

#Create data directory for treatment data
data_t.dir <- "data path address"
raw_t_data <- Read10X(data.dir=data_t.dir)
list.files(data_t.dir)

#Create Seurat object for treatment 
seurat_obj2 <-CreateSeuratObject(
  counts = raw_t_data,
  project = 'Arabidopsis_t',
  min.cells = 3, # genes detected in >3 cells
  min.features = 300 #cell with >200 genes
)
seurat_obj2
head(colnames(seurat_obj2))
#Add metadata
seurat_obj2$sample <-'treatment'
seurat_obj2$condition <- 'treatment'
head(seurat_obj2@meta.data)

#Add mitochondrial Percentage data
seurat_obj2[["percent.mt"]] <- PercentageFeatureSet(
  seurat_obj2,
  pattern = "^ATMG"
)
#Add Chlorplast Percentage data
seurat_obj2[["percent.chloro"]] <- PercentageFeatureSet(
  seurat_obj2,
  pattern = "^ATCG"
)

#Merge both the object
seurat_obj <- merge(
  seurat_obj1,
  y = seurat_obj2,
  add.cell.ids = c("CTRL", "TRT")
)

#Filter the low quality data
#Cell expressing less than 200 genes (likely dead), more than 6000 genes(potential doublets), 
#>10% mitochondria gene content(stress/dead)

seurat_obj <-subset(seurat_obj,
                    subset =nFeature_RNA > 300 &
                      nFeature_RNA < 6000 &
                      percent.mt < 10 &
                      percent.chloro < 20)

#To visualize key QC matics per cell using violin plots:
#-nFeatures_RNA(genes/cell), ncount_RNA(rna/cell), percent.mt, percent.chloro
VlnPlot(seurat_obj,
        features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt', 'percent.chloro'),
        pt.size = 0.1,
        ncol = 3)

#Create the scatter plot
FeatureScatter(
  seurat_obj,
  feature1="nCount_RNA",
  feature2="nFeature_RNA"
)
head(colnames(seurat_obj))

#log-Normalization, transform raw UMI into normalised values to account for sequencing depth
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

#Identify 2000 highly variable genes
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = 'vst', nfeatures = 2000)

#Visualize the variable features
var_plot <- VariableFeaturePlot(seurat_obj)
LabelPoints(plot = var_plot, points = head(VariableFeatures(seurat_obj), 20), repel = TRUE)
VariableFeaturePlot(seurat_obj)
length(VariableFeature(seurat_obj))

#Dimentionality reduction & clustering, scaling centers gene expression
seurat_obj <-ScaleData(seurat_obj)

#Run PCA to reduce dim and identify major source of variation
#npcs = number of Principal components to compute (here we are using 50)
seurat_obj <-RunPCA(seurat_obj, npcs = 50)

#Visualize the SD of each PC to decide how many PCs to retain
#A sharp drop("elbow") in the plot helps determine the optimal number of PCs
ElbowPlot(seurat_obj, ndims=50)

#Based on the elbow plot, we select first 13 PC 
pcs <- 13

#Identify the k-nearest neighbor for each cell based on the selected PCs
seurat_obj <-FindNeighbors(seurat_obj, dims = 1:pcs)

#Perform graph based culturing of the cells
#Resolution affects the granularity: higher = more clusters
seurat_obj <-FindClusters(seurat_obj, resolution = 0.9)

#Run UMAP for 2D visualization of the clusters
seurat_obj <-RunUMAP(seurat_obj, dims = 1:pcs)

#Plot UMAP clustering result
DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE) + 
  ggtitle("Sc-profilling of Arabidopsis leaf")

#Cell type annotation
#Manually create .csv reference file 
library(readr)
marker_ref <- read_csv("Arabidopsis_leaf_marker.csv", show_col_types = FALSE)
head(colnames(marker_ref))

#Join the layer in seurat_obj
Layers(seurat_obj)
seurat_obj <- JoinLayers(seurat_obj)
Layers(seurat_obj)

#Find marker genes for all clusters
#Only return genes with higher expression in the cluster compared to others
#only.pos= only positive, min.pct=expressed in atleast 25% of cell in a cluster
#min log-fold change
markers <- FindAllMarkers(seurat_obj,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)
write.csv(markers, "markers.csv", row.names = FALSE)

#Match your marker genes with reference genes
matched_markers <- merge(
  markers,
  marker_ref,
  by.x="gene",
  by.y="GeneID"
)

#Count marker matches per cluster
library(dplyr)
annotation_table <- matched_markers %>%
  group_by(cluster, Gene_type) %>%
  summarise(
    n_markers=n()
  ) %>%
  arrange(cluster, desc(n_markers))

#Add annotation to seurat
cluster_identity <- c(
  "0"="Bundle sheath/Hydathode",
  "1"="Mesophyll/Bundle sheath",
  "2"="Mesophyll",
  "3"="Bundle sheath/Hydathode",
  "4"="Mesophyll",
  "5"="Mesophyll",
  "6"="Bundle sheath/Hydathode",
  "7"="Mesophyll",
  "8"="Mesophyll/Bundle sheath",
  "9"="Mesophyll",
  "10"="Mesophyll",
  "11"="Mesophyll/Bundle sheath/Parenchyma",
  "12"="Mesophyll",
  "13"="Mesophyll/Bundle sheath",
  "14"="Bundle sheath/Hydathode",
  "15"="Companion cell/Hydathode",
  "16"="Parenchyma/Hydathode",
  "17"="Mesophyll",
  "18"="Gaurd cell/Procambium"
)
levels(Idents(seurat_obj))
colnames(seurat_obj@meta.data)
Idents(seurat_obj) <- "seurat_clusters"

#Rename clusters
seurat_obj <- RenameIdents(
  seurat_obj,
  cluster_identity
)
#Save the annotations in the metadata
seurat_obj$Gene_type <- Idents(seurat_obj)

#Plot cell type in UMAP
DimPlot(seurat_obj, group.by = "Gene_type", label = TRUE)

#Identify top markers
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(n =5, order_by = avg_log2FC)

write.csv(top_markers, "highexp_markers.csv")

#Current indents are cell annotation change to cluster
table(seurat_obj$seurat_clusters)
Idents(seurat_obj) <- "seurat_clusters"
levels(Idents(seurat_obj))

#Differential expression between two specific cell type
de_genes <-FindMarkers(seurat_obj,
                       ident.1 = 4,
                       ident.2 = 14,
                       min.pct =0.25,
                       group.by = "seurat_clusters",
                       logfc.threshold =0.25)
#View the top DE genes
head(de_genes, 8)
#Save to csv
write.csv(de_genes, "de_cluster4vs14.csv")

#Add row for significance
de_genes$gene <-rownames(de_genes)
de_genes$significant <- ifelse(de_genes$p_val_adj <0.05 & abs(de_genes$avg_log2FC) >0.5, "Yes", "No")

#Plot using ggplot2
library(ggplot2)
ggplot(de_genes, aes(x= avg_log2FC, y= -log10(p_val_adj), 
                         colour = significant)) +
  geom_point(alpha= 0.8) +
  scale_color_manual(values = c('grey' , 'darkgreen')) +
  theme_minimal() +
  labs(title = "Volcano plot: cluster4 vs cluster14",
       x= 'Log2 fold Change',
       y= "-Log10 Adjusted P-value")

#Differential expression of gene in between control vs treatment
Idents(seurat_obj) <- "condition"
de_genes1 <-FindMarkers(seurat_obj,
                        ident.1 = "treatment",
                        ident.2 = "control",
                        min.pct =0.25,
                       logfc.threshold = 0.25)
de_genes1$gene <-rownames(de_genes1)
de_genes1$significant <- ifelse(de_genes1$p_val_adj <0.05 & abs(de_genes1$avg_log2FC) >0.5, "Yes", "No")

#Annotate which cluster do each gene belong to
library(dplyr)
gene_clusters <- markers %>%
  group_by(gene) %>%
  summarise(
    cluster = paste(unique(cluster), collapse = ", "),
    .groups = "drop"
  )

de_genes1 <- de_genes1 %>%
  mutate(gene = rownames(de_genes1)) %>%
  left_join(gene_clusters, by = "gene")

#View the top DE genes
head(de_genes1, 8)
#save to csv
write.csv(de_genes1, "de_treatmentvscontrol.csv")

#Since most of the p_val is 0 , it is replaced with small value
de_genes$p_val_adj_plot <- pmax(de_genes$p_val_adj, 1e-50)

#Plot using ggplot2
library(ggplot2)
ggplot(de_genes1, aes(x= avg_log2FC, y= -log10(p_val_adj_plot), 
                     colour = significant)) +
  geom_point(alpha= 0.8) +
  scale_color_manual(values = c('grey' , 'darkblue')) +
  theme_minimal() +
  labs(title = "Volcano plot: treatment vs control",
       x= 'Log2 fold Change',
       y= "-Log10 Adjusted P-value")


#Differential expression of gene in between control vs treatment in Mesophyll cells

#Check the diff expression based on the cell type
Idents(seurat_obj) <- "Gene_type"
Mesophyll_diff <- subset(seurat_obj, idents = "Mesophyll")

Idents(Mesophyll_diff) <- "condition"

de_mesophyll <- FindMarkers(
  Mesophyll_diff,
  ident.1 = "treatment",
  ident.2 = "control"
)

#View the top DE genes
head(de_mesophyll, 8)
#Save to csv
write.csv(de_mesophyll, "de_mesophyll_controlvs_treatment.csv")

#Add column for significance
de_mesophyll$gene <-rownames(de_mesophyll)
de_mesophyll$significant <- ifelse(de_mesophyll$p_val_adj <0.05 & abs(de_mesophyll$avg_log2FC) >0.5, "Yes", "No")

#Plot using ggplot2
library(ggplot2)
ggplot(de_mesophyll, aes(x= avg_log2FC, y= -log10(p_val_adj), 
                      colour = significant)) +
  geom_point(alpha= 0.8) +
  scale_color_manual(values = c('grey' , 'red')) +
  theme_minimal() +
  labs(title = "Volcano plot: control vs treatment",
       x= 'Log2 fold Change',
       y= "-Log10 Adjusted P-value")

