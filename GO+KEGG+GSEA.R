#GO+KEGG+GSEA
library(openxlsx)
library(ggplot2)
library(stringr)
library(enrichplot)
library(clusterProfiler)
library(GOplot)
library(DOSE)
library(ggnewscale)
library(topGO)
library(circlize)
library(ComplexHeatmap)



deg<-sig_dge.cluster

colnames(deg)
logFC=0
P.Value = 0.05
type1 = (deg$p_val< P.Value)&(deg$avg_log2FC< -logFC)
type2 = (deg$p_val< P.Value)&(deg$avg_log2FC> -logFC)
deg$Group2 = ifelse(type1,"Down",ifelse(type2,"Up","Not-Sig"))
table(deg$Group2)


diff.genes<-rownames(subset(deg,Group2!='Not-Sig'))


diff.df <- bitr(diff.genes,
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)

gene<-diff.df
GO_database <- 'org.Hs.eg.db'
GO<-enrichGO( gene$ENTREZID,
              OrgDb = GO_database,
              keyType = "ENTREZID",
              ont = "ALL",
              pvalueCutoff = 1,
              qvalueCutoff = 1,
              readable = T)
GO

R.utils::setOption("clusterProfiler.download.method",'auto')

KEGG_database <- 'hsa'
KEGG<-enrichKEGG(gene$ENTREZID,
                 organism = KEGG_database,
                 pvalueCutoff = 1,
                 qvalueCutoff = 1)

KEGG

names(deg) <- c('SYMBOL','Log2FoldChange')
info_merge <- merge(deg,gene,by='SYMBOL')
GSEA_input <- info_merge$Log2FoldChange
names(GSEA_input) = info_merge$ENTREZID
GSEA_input = sort(GSEA_input, decreasing = TRUE)
GSEA_KEGG <- gseKEGG(GSEA_input, organism = KEGG_database, pvalueCutoff = 1)#GSEA富集分析
enrichplot::cnetplot(GO,circular=FALSE,colorEdge = TRUE)
enrichplot::cnetplot(KEGG,circular=FALSE,colorEdge = TRUE)

outFile="Rplot.pdf"
pdf(file=outFile,width=12,height=8)
enrichplot::cnetplot(GO,circular=FALSE,colorEdge = TRUE)
dev.off()











ego_ALL <- enrichGO(gene          = row.names(sig_dge.cluster),
                    #universe     = row.names(dge.celltype),
                    OrgDb         = 'org.Hs.eg.db',
                    keyType       = 'SYMBOL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1)
ego_all <- data.frame(ego_ALL)
write.csv(ego_ALL,'AS_neipi_enrichGO.csv')

ego_CC <- enrichGO(gene          = row.names(sig_dge.cluster),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)

ego_MF <- enrichGO(gene          = row.names(sig_dge.cluster),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)

ego_BP <- enrichGO(gene          = row.names(sig_dge.cluster),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)

ego_CC@result$Description <- substring(ego_CC@result$Description,1,70)
ego_MF@result$Description <- substring(ego_MF@result$Description,1,70)
ego_BP@result$Description <- substring(ego_BP@result$Description,1,70)

p_BP <- barplot(ego_BP,showCategory = 10) + ggtitle("barplot for Biological process")
p_CC <- barplot(ego_CC,showCategory = 10) + ggtitle("barplot for Cellular component")
p_MF <- barplot(ego_MF,showCategory = 10) + ggtitle("barplot for Molecular function")
plotc <- p_BP/p_CC/p_MF


remotes::install_github("YuLab-SMU/createKEGGdb")
createKEGGdb::create_kegg_db('hsa')
install.packages("KEGG.db_1.0.tar.gz",repos=NULL, type = "source")
library(KEGG.db)
library(org.Hs.eg.db)
library(clusterProfiler)
genename <- rownames(sig_dge.cluster)
View(genename)
gene <-bitr(genename, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db",drop= TRUE)
enrich.KEEG <- enrichKEGG( gene = gene$ENTREZID,
                           organism = "hsa",
                           keyType = "kegg",
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH",
                           use_internal_data = FALSE)
keeg <- as.data.frame(KEGG_enrich)
barplot(keeg)
dotplot(keeg)
browseKEGG(keeg, 'hsa04110')

ggplot(KEGG_enrich)










