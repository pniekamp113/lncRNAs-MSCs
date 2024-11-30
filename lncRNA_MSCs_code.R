#Understanding the role of lncRNAs in MSC stemness and differentiation


#Load required libraries 


library(dplyr)
library(Seurat)
library(patchwork)
library(readr)
library(RColorBrewer)
library(ggplot2)
library(VisCello)
library(viridis)
library(forcats)
library(ggrastr)
library(cowplot)
library(GSEABase)
library(SeuratDisk)
library(ComplexHeatmap)
library(irlba)
library(AUCell)
library(biomaRt)
library(glmGamPoi)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(glmGamPoi)


#load seurat object
getwd()
MSC_final <- readRDS("C:\\Users\\Patrick Niekamp\\Desktop\\Experiments\\Exp. 57_Cell paper scRNA seq\\MSC_final.rds")

DimPlot(MSC_final, label = TRUE)



###########
#check for lncRNAs

head(rownames(MSC_final)) #Gene symbols

# Connect to Ensembl database
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl") 

# Retrieve gene annotations
gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "gene_biotype", "external_gene_name"),
  filters = "external_gene_name",
  values = rownames(MSC_final),
  mart = mart
)

head(gene_mapping)
unique(gene_mapping$gene_biotype)


# Filter for lncRNAs
lncRNA_genes <- gene_mapping[gene_mapping$gene_biotype == "lncRNA", ]
head(lncRNA_genes)

#Filter seurat object to focus on lncRNAs
lncRNA_features <- intersect(lncRNA_genes$external_gene_name, rownames(MSC_final))
print(lncRNA_features)

MSC_lncRNA <- subset(MSC_final, features = lncRNA_features)
saveRDS(MSC_lncRNA, file = "MSC_final_lncRNA.rds")


lncRNA_markers <- FindAllMarkers(MSC_lncRNA, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)
head(lncRNA_markers)
write.csv(lncRNA_markers, file = "lncRNA_MSC_markers.csv")


#show expression of top 15 genes per cluster using dotplot

top_10 <- lncRNA_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 10)

top_10 <- top_10$gene
top_10 <- unique(top_10)

DotPlot(MSC_lncRNA, features = top_10) +
  scale_color_gradientn(colors = c("grey", "blue")) +  # Adjust color gradient as needed
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10)) +
  ggtitle("MSC markers")


#####
# Use 