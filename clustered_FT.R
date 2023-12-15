library(Seurat)
library(patchwork)
library(cowplot)
library(data.table)
library(SeuratWrappers)
library(SeuratData)
library(harmony)
library(dplyr)

#1.0 resolution clustering
FTC_GEX1 <- readRDS('/home/jurny0402/FTC/code/saved_codes/reclustering/230502_FTC_GEX_reclustered_res1.0.rds')
View(FTC_GEX1@meta.data)
DimPlot(FTC_GEX1, label = TRUE, label.size = 4)
DimPlot(FTC_GEX1, label = TRUE, label.size = 4, split.by = "tissue_type")
#reclustering
FTC_GEX1_final <- WhichCells(FTC_GEX1, ident = c(0:3,5:16,18:23,25:27,29:30,36,40,42))
FTC_GEX2 <- subset(FTC_GEX1, cells=FTC_GEX1_final)

DimPlot(FTC_GEX2, label = T, label.size = 4, split.by = "tissue_type")
DimPlot(FTC_GEX2, group.by = 'tissue_type')
#FT cell들만 골라내기
FT <- subset(FTC_GEX2, tissue_type == "tumor")
DimPlot(FT, label = T, label.size = 4)
DimPlot(FT, label = T, label.size = 4, split.by = "tissue_type")

saveRDS(FT, file = '/home/jurny0402/TCR_explore/code/clustered_FT.rds')
save(FT, file = '/home/jurny0402/TCR_explore/code/clustered_FT.rda')
table(FT$seurat_clusters)

view <- readRDS("/home/jurny0402/TCR_explore/code/contigid_barcode_adjusted_clustered_FT.rds")
View(view@meta.data)
rownames(view@meta.data)[1]
----------------------------------------------------------
FTC_GEX1$orig.ident
adenoma <- subset(FTC_GEX2, subset = orig.ident == "FT1"|orig.ident == "FT2"|orig.ident == "FT3")
carcinoma <- subset(FTC_GEX2, subset = orig.ident == "FT4"|orig.ident == "FT5")

DimPlot(adenoma, label = T, label.size = 4)
DimPlot(carcinoma, label = T, label.size = 4)

saveRDS(adenoma, file = '/home/jurny0402/TCR_explore/code/230530_adenoma.rds')
saveRDS(carcinoma, file = '/home/jurny0402/TCR_explore/code/230530_carcinoma.rds')


FT1 <- read.csv("/node200data/FTC/data/cellranger_multi_outs/FT1/outs/per_sample_outs/FT1/vdj_t/filtered_contig_annotations.csv")
FT2 <- read.csv("/node200data/FTC/data/cellranger_multi_outs/FT2/outs/per_sample_outs/FT2/vdj_t/filtered_contig_annotations.csv")
FT3 <- read.csv("/node200data/FTC/data/cellranger_multi_outs/FT3/outs/per_sample_outs/FT3/vdj_t/filtered_contig_annotations.csv")
FT4 <- read.csv("/node200data/FTC/data/cellranger_multi_outs/FT4/outs/per_sample_outs/FT4/vdj_t/filtered_contig_annotations.csv")
FT5 <- read.csv("/node200data/FTC/data/cellranger_multi_outs/FT5/outs/per_sample_outs/FT5/vdj_t/filtered_contig_annotations.csv")

View(FT1)

