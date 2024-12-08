library(monocle)
packageVersion("igraph")
packageVersion('monocle')
Hu_AO_db_QC2<-sub_pbmc02
Hu_AO_db_QC2$New_celltype<-Hu_AO_db_QC2$celltype

Idents(Hu_AO_db_QC2) <-Hu_AO_db_QC2$New_celltype
DimPlot(Hu_AO_db_QC2, reduction = "umap",group.by = "New_celltype",label = T)

SMC_cells <- subset(Hu_AO_db_QC2, idents = c("SMC1", "SMC2", "SMC3", "SMC4"))
DimPlot(SMC_cells,reduction = "umap",group.by = "New_celltype",label = T)
SMC_cells<-Hu_AO_db_QC2

monocle.matrix <- as(GetAssayData(object = SMC_cells, layer = "counts", assay = "RNA"), 'sparseMatrix')
monocle.sample <- as.matrix(SMC_cells@meta.data)
monocle.geneAnn <- data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
monocle.sample <- as.data.frame(monocle.sample)
monocle.geneAnn <- as.data.frame(monocle.geneAnn)


pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- newCellDataSet(monocle.matrix, phenoData = pd, featureData = fd)


cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1)
cds <- cds[fData(cds)$num_cells_expressed > 10, ]


expressed_genes <- rownames(subset(fData(cds), num_cells_expressed >= 10))

diff_test_res <- differentialGeneTest(cds[expressed_genes, ], fullModelFormulaStr = "~New_celltype")

ordering_genes <- rownames(subset(diff_test_res, qval < 0.01))

cds <- setOrderingFilter(cds, ordering_genes)


cds <- reduceDimension(cds, method = 'DDRTree')

cds <- orderCells(cds)
plot_ordering_genes(cds)
plot_cell_trajectory(cds, color_by = "New_celltype")
plot_cell_trajectory(cds, color_by = "Pseudotime")
plot_cell_trajectory(cds, color_by = "State")
plot_cell_trajectory(cds, color_by = "State")+facet_wrap(~State, nrow = 3)
plot_cell_trajectory(cds, color_by = "State") +
  facet_wrap("~State", nrow = 3)
plot_cell_trajectory(cds, color_by = "New_celltype") +
  facet_wrap("~New_celltype", nrow = 3)
plot_cell_trajectory(cds, color_by = "orig.ident") +
  facet_wrap("~orig.ident", nrow = 3)
plot_cell_trajectory(cds, color_by = "orig.ident") +
  facet_wrap("~orig.ident", nrow = 3)+ scale_color_manual(values = c("#d86967", "#58539f"))
plot_complex_cell_trajectory(cds,color_by = "New_celltype") +  theme(legend.title = element_blank())  # 去除图例名
cg<-as.character(ordering_genes[1:3])
plot_genes_in_pseudotime(cds[cg,],color_by = "New_celltype")
p <- plot_genes_jitter(cds[cg,],                       
                       grouping = "New_celltype",                       
                       color_by = "New_celltype",
                       nrow=3,                       
                       ncol = NULL)
p + theme(axis.text.x = element_blank())
library(ggsci)
pData(cds)$ITGAM=log2(exprs(cds)['ITGAM',]+1)
plot_cell_trajectory(cds, color_by = 'ITGAM') +
  scale_color_gsea()

library(ggsci)
pData(cds)$GP1BA=log2(exprs(cds)['GP1BA',]+1)
plot_cell_trajectory(cds, color_by = 'GP1BA') +
  scale_color_gsea()

plot_cell_trajectory(cds, color_by = 'GP1BA') +
  scale_color_gradient(low = 'grey', high = 'red')

library(patchwork)
topN <- head(ordering_genes, 3)
topN <- c("ITGAM", "GP1BA")
topN_cds <- cds[topN,]
p1 <-plot_genes_in_pseudotime(topN_cds, color_by = "Pseudotime")
p2 <-plot_genes_in_pseudotime(topN_cds, color_by = "celltype")