library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)


dir = c('pso1/', "pso2/")
names(dir) = c('PSO', "PSO")
counts <- Read10X(data.dir =dir)
scRNA=CreateSeuratObject(counts,min.cells = 3,project="os",min.features = 300)

scRNA[["mt_percent"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
HB_genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB_genes, rownames(scRNA@assays$RNA))
HB_genes <- rownames(scRNA@assays$RNA)[HB_m]
HB_genes <- HB_genes[!is.na(HB_genes)]
scRNA[["HB_percent"]] <- PercentageFeatureSet(scRNA, features=HB_genes)

PSO <- subset(scRNA,
              subset = nFeature_RNA > 300 & nFeature_RNA < 7000 &
                mt_percent < 10 &
                HB_percent < 3 &
                nCount_RNA < quantile(nCount_RNA,0.97) & nCount_RNA > 1000)


dir = c('aspso1/','aspso2/')
names(dir) = c('ASPSO','ASPSO')
counts <- Read10X(data.dir =dir)
scRNA=CreateSeuratObject(counts,min.cells = 3,project="os",min.features = 300)

scRNA[["mt_percent"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
HB_genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB_genes, rownames(scRNA@assays$RNA))
HB_genes <- rownames(scRNA@assays$RNA)[HB_m]
HB_genes <- HB_genes[!is.na(HB_genes)]
scRNA[["HB_percent"]] <- PercentageFeatureSet(scRNA, features=HB_genes)

Aspso <- subset(scRNA,
                subset = nFeature_RNA > 300 & nFeature_RNA < 7000 &
                  mt_percent < 10 &
                  HB_percent < 3 &
                  nCount_RNA < quantile(nCount_RNA,0.97) & nCount_RNA > 1000)

Aspso$batch <- 'ASPSO'
PSO$batch <- 'PSO'


library(harmony)
combined <- merge(Aspso, y = list(PSO), add.cell.ids = c("ASPSO", "PSO"))
#combined<-Pso
combined <- NormalizeData(combined, verbose = FALSE)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, features = VariableFeatures(object = combined))
head(combined@meta.data)

harmony_integrated <- RunHarmony(combined, group.by.vars = "batch", plot_convergence = TRUE)
scRNA1<-harmony_integrated
#scRNA1<-combined
pc.num=1:5
scRNA1<-FindNeighbors(scRNA1,dims = pc.num)
scRNA1<-FindClusters(scRNA1,resolution = 0.8)
scRNA1<-BuildClusterTree(scRNA1)

# UMAP
scRNA1 <- RunUMAP(scRNA1, reduction = "harmony", dims = 1:5)


library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggunchull)
library(tidydr)
library(ggsci)
library(Cairo)

color1 <- c("#FFBE7A","#FA7F6F",
            "#63b2ee","#76da91","#7CAE00",
            "#f89588","#9192ab","#C77CFF","#FFA500")

DimPlot(scRNA1, reduction = "umap",#group.by = "orig.ident",
        #cols = color1,
        pt.size = 0.5,
        label = F,label.box = F
)+theme_dr(xlength = 0.22, ylength = 0.22,
           arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())


celltypemakers<-c("CDH5","ICAM1","PECAM1",
                  "COL1A1","COL1A2","COL3A1",
                  "HLA-DRB1","HLA-DQA1","CD68",
                  "CD2","CD3D","IL32")

DotPlot(scRNA1, features = celltypemakers)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))


cluster_types <- data.frame(
  cluster = c(0,1,2,3,4,5,6,7,8,9,10,
              11,12,13,14,15,16,17,18,19,20,21),
  celltype = c("T_cells", "Fibroblasts", "Endothelial_cells","Endothelial_cells","Eccrine_gland","Mast_cells",
               "T_cells","Fibroblasts","Fibroblasts","Fibroblasts","Endothelial_cells",
               "Keratinocytes","Macrophages","T_cells","Endothelial_cells","Keratinocytes",
               "DC","Macrophages","Keratinocytes","Fibroblasts","Endothelial_cells",
               "Endothelial_cells")
)

scRNA1@meta.data$celltype ="NA"
for(i in 1:nrow(cluster_types)) {
  cluster_id <- cluster_types$cluster[i]
  cell_type <- cluster_types$celltype[i]
  scRNA1$celltype[scRNA1$seurat_clusters == cluster_id] <- cell_type
}

cluster_types <- data.frame(
  cluster = c(0,1,2,3,4,5,6,7,8,9,10,
              11,12,13,14,15,16,17,18,19,20,21),
  celltype = c("T_cells", "Fibroblasts", "Endothelial_cells","Endothelial_cells","Eccrine_gland","Mast_cells",
               "T_cells","Fibroblasts","Fibroblasts","Fibroblasts","Endothelial_cells",
               "Keratinocytes","Monocyte-derived_Macrophages","T_cells","Endothelial_cells","Keratinocytes",
               "DC","Tissue-resident_Macrophages","Keratinocytes","Fibroblasts","Endothelial_cells",
               "Endothelial_cells")
)

scRNA1@meta.data$celltype ="NA"
for(i in 1:nrow(cluster_types)) {
  cluster_id <- cluster_types$cluster[i]
  cell_type <- cluster_types$celltype[i]
  scRNA1$celltype[scRNA1$seurat_clusters == cluster_id] <- cell_type
}

DimPlot(scRNA1, reduction = "umap",#group.by = "orig.ident",
        #cols = color1,
        pt.size = 0.5,
        label = T,label.box = F
)+theme_dr(xlength = 0.22, ylength = 0.22,
           arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())

DimPlot(scRNA1, reduction = "umap",group.by = "celltype",
        #cols = color1,
        pt.size = 0.5,
        label = F,label.box = F
)+theme_dr(xlength = 0.22, ylength = 0.22,
           arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())








pbmc1 <- scRNA1
sub_pbmc02_matrix<- pbmc1@assays$RNA@counts[,which(pbmc1@meta.data$seurat_clusters%in%c(12,17))]


sub_pbmc02 <- CreateSeuratObject(counts = sub_pbmc02_matrix,project = 'sub_02')
sub_pbmc02[['percent.mt']] <- PercentageFeatureSet(sub_pbmc02,pattern = '^MT-')
sub_pbmc02 <- NormalizeData(sub_pbmc02,normalization.method = "LogNormalize",scale.factor = 10000)
sub_pbmc02 <- FindVariableFeatures(sub_pbmc02,selection.method = 'vst',nfeatures = 2000)
sub_pbmc02 <- ScaleData(sub_pbmc02,vars.to.regress = "percent.mt")
sub_pbmc02<-RunPCA(sub_pbmc02,features = VariableFeatures(scRNA1))
DimPlot(sub_pbmc02,reduction = "pca",group.by = "orig.ident")
ElbowPlot(sub_pbmc02,ndims=5,reduction = "pca")

sub_pbmc02 <- RunHarmony(sub_pbmc02, group.by.vars = "orig.ident", plot_convergence = TRUE)

sub_pbmc02 <- FindNeighbors(sub_pbmc02, dims = 1:5)
sub_pbmc02 <- FindClusters(sub_pbmc02, resolution = 1)
head(Idents(sub_pbmc02),8)
sub_pbmc02<- RunUMAP(sub_pbmc02, dims = 1:10)


cluster_types <- data.frame(
  cluster = c(0,1,2,3,4,5,6,7,8,9,10),
  celltype = c("Tissue-resident_Macrophages", "Monocyte-derived_Macrophages", "Monocyte-derived_Macrophages","Tissue-resident_Macrophages","Monocyte-derived_Macrophages","Tissue-resident_Macrophages",
               "Tissue-resident_Macrophages","Monocyte-derived_Macrophages","Monocyte-derived_Macrophages","Tissue-resident_Macrophages","Monocyte-derived_Macrophages")
)

sub_pbmc02@meta.data$celltype ="NA"
for(i in 1:nrow(cluster_types)) {
  cluster_id <- cluster_types$cluster[i]
  cell_type <- cluster_types$celltype[i]
  sub_pbmc02$celltype[sub_pbmc02$seurat_clusters == cluster_id] <- cell_type
}



color1<-c("#00BFC4","#F8766D")

DimPlot(sub_pbmc02, reduction = "umap",group.by = "orig.ident",
        #cols = color1,
        pt.size = 0.5,
        label = T,label.box = F
)+theme_dr(xlength = 0.22, ylength = 0.22,
           arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())

DimPlot(sub_pbmc02, reduction = "umap",group.by = "celltype",
        #cols = color1,
        pt.size = 0.5,
        label = F,label.box = F
)+theme_dr(xlength = 0.22, ylength = 0.22,
           arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())









