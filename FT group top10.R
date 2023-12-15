library(Seurat)
library(patchwork)
library(cowplot)
library(data.table)
library(SeuratWrappers)
library(SeuratData)
library(harmony)
library(dplyr)
library(SingleR)

#Normalization된 자료
FTC_GEX <- readRDS('/home/jurny0402/FTC/code/saved_codes/230404_FTC_GEX_SCTranformed_(beforeFastMNN).rds')


Seurat <- readRDS("/home/jurny0402/TCR_explore/code/contigid_barcode_adjusted_clustered_FT.rds")
View(Seurat@meta.data)
View(FTC_GEX1@meta.data)
