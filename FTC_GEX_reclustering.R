library(Seurat)
library(patchwork)
library(cowplot)
library(data.table)
library(SeuratWrappers)
library(SeuratData)
library(harmony)
library(dplyr)
library(SingleR)
  devtools::install_github('dviraran/SingleR')
  #install.packages("remotes")
  #remotes::install_github('dviraran/SingleR')
  
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("SingleR")
  
library(repr)
library(batchelor)
library(glmGamPoi)
library(stringr)
library(tibble)
library(DoubletFinder)
  install.packages("remotes")
  remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
library(ggplot2)

file_list<-list.files("/home/jurny0402/FTC/code/221124_GEX_count_matrices_after_cellbender_scrublet")
file_list
file_list <- file_list[-3]
file_list

meta.info=data.table(file_path = file_list)
meta.info$sampleID <- c("FT1", 'FT2', 'FT3', 'FT4','FT5', 'Thy1', 'Thy10', 'Thy4', 'Thy5', 'Thy6', 'Thy7')
DT::datatable(meta.info)
meta.info$file_path
meta.info$sampleID[i]

sample_paths = c("/home/jurny0402/FTC/code/221124_GEX_count_matrices_after_cellbender_scrublet/FT1_count_matrix.csv",
                 "/home/jurny0402/FTC/code/221124_GEX_count_matrices_after_cellbender_scrublet/FT2_count_matrix.csv",
                 "/home/jurny0402/FTC/code/221124_GEX_count_matrices_after_cellbender_scrublet/FT3_count_matrix.csv",
                 "/home/jurny0402/FTC/code/221124_GEX_count_matrices_after_cellbender_scrublet/FT4_count_matrix.csv",
                 "/home/jurny0402/FTC/code/221124_GEX_count_matrices_after_cellbender_scrublet/FT5_count_matrix.csv",
                 "/home/jurny0402/FTC/code/221124_GEX_count_matrices_after_cellbender_scrublet/Thy1_count_matrix.csv",
                 "/home/jurny0402/FTC/code/221124_GEX_count_matrices_after_cellbender_scrublet/Thy10_count_matrix.csv",
                 "/home/jurny0402/FTC/code/221124_GEX_count_matrices_after_cellbender_scrublet/Thy4_count_matrix.csv",
                 "/home/jurny0402/FTC/code/221124_GEX_count_matrices_after_cellbender_scrublet/Thy5_count_matrix.csv",
                 "/home/jurny0402/FTC/code/221124_GEX_count_matrices_after_cellbender_scrublet/Thy6_count_matrix.csv",
                 "/home/jurny0402/FTC/code/221124_GEX_count_matrices_after_cellbender_scrublet/Thy7_count_matrix.csv")

batch_list = list()

for(i in 1:5) {
  dir_of_10X = meta.info$file_path[i]
  sampleID = meta.info$sampleID[i]
  meta_mat = read.csv(sample_paths[i], row.names = 1)
  FTC_GEX <- CreateSeuratObject(counts = meta_mat, project = sampleID, min.cells = 3, min.features = 200)
  FTC_GEX[["percent.mt"]] <- PercentageFeatureSet(FTC_GEX, pattern = "^MT-")
  FTC_GEX@meta.data[,'SampleID'] = sampleID
  FTC_GEX@meta.data[,'tissue_type'] = 'tumor'
  batch_list <- append(batch_list, FTC_GEX)
  print(paste("Finishing processing", sampleID, 'at', Sys.time()))
}

for(i in 6:nrow(meta.info)) {
  dir_of_10X = meta.info$file_path[i]
  sampleID = meta.info$sampleID[i]
  meta_mat = read.csv(sample_paths[i], row.names = 1)
  FTC_GEX <- CreateSeuratObject(counts = meta_mat, project = sampleID, min.cells = 3, min.features = 200)
  FTC_GEX[["percent.mt"]] <- PercentageFeatureSet(FTC_GEX, pattern = "^MT-")
  FTC_GEX@meta.data[,'SampleID'] = sampleID
  FTC_GEX@meta.data[,'tissue_type'] = 'normal'
  batch_list <- append(batch_list, FTC_GEX)
  print(paste("Finishing processing", sampleID, 'at', Sys.time()))
}

FTC_GEX <- Reduce(function(x,y){merge(x,y,merge.data = TRUE)}, batch_list)

saveRDS(FTC_GEX, file = '/home/jurny0402/FTC/code/saved_codes/230404_FTC_GEX_merged_(beforeqc_noSCT).rds')
FTC_GEX <- readRDS('/home/jurny0402/FTC/code/saved_codes/230404_FTC_GEX_merged_(beforeqc_noSCT).rds')
#맨 처음 시작시, workspace loading 후 위와 같이 저장된 코드 가져오기

VlnPlot(FTC_GEX, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, combine = F)
VlnPlot(FTC_GEX, features = "percent.mt", pt.size = 0.1, y.max = 20)

FTC_GEX <- subset(FTC_GEX, subset = percent.mt <= 10 & nFeature_RNA >= 500 & nFeature_RNA <= 10000 & nCount_RNA >= 1000 & nCount_RNA <= 10000)
VlnPlot(FTC_GEX, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)

table(FTC_GEX@meta.data$SampleID)
table(FTC_GEX@meta.data$tissue_type)

FTC_GEX <- SCTransform(FTC_GEX, vars.to.regress = c('nCount_RNA', 'percent.mt'), vst.flavor = 'v2', method = 'glmGamPoi', verbose = F)
saveRDS(FTC_GEX, file = '/home/jurny0402/FTC/code/saved_codes/230404_FTC_GEX_SCTranformed_(beforePCA).rds')
FTC_GEX <- readRDS('/home/jurny0402/FTC/code/saved_codes/230404_FTC_GEX_SCTranformed_(beforeFastMNN).rds')
#중간저장

FTC_GEX <- RunPCA(FTC_GEX)
saveRDS(FTC_GEX, file = '/home/jurny0402/FTC/code/saved_codes/230407_FTC_GEX_afterpca_(beforeFastMNN)')
FTC_GEX <- readRDS('/home/jurny0402/FTC/code/saved_codes/230407_FTC_GEX_afterpca_(beforeFastMNN)')
#중간저장

ElbowPlot(FTC_GEX, ndims=50)
DimHeatmap(FTC_GEX, dims=1:15,cells=5000)

#harmony 안 할 거라서 FTC_GEX1 지정 안하고 할라했는데 DietSeurat 필요
#FTC_GEX <- RunFastMNN(object.list = SplitObject(object = FTC_GEX, split.by = "SampleID"), verbose = F)

FTC_GEX1 <- DietSeurat(FTC_GEX)
FTC_GEX1 <- RunFastMNN(object.list = SplitObject(FTC_GEX1, split.by = "SampleID"), verbose = F)

saveRDS(FTC_GEX1, file = '/home/jurny0402/FTC/code/saved_codes/230407_FTC_GEX_fastMNN_(stringentqc_SCT).rds')
FTC_GEX1 <- readRDS('/home/jurny0402/FTC/code/saved_codes/230407_FTC_GEX_fastMNN_(stringentqc_SCT).rds')
#중간저장

DimPlot(FTC_GEX, reduction="pca")
DimPlot(FTC_GEX1, reduction="mnn")

num_pc = 40
#pca-UMAP(dims=1:40)
FTC_GEX <- RunUMAP(FTC_GEX, reduction = 'pca', dims = 1:num_pc, verbose = F)
DimPlot(FTC_GEX, group.by = 'SampleID')
DimPlot(FTC_GEX, group.by = 'tissue_type')

#mnn-UMAP(dims=1:40)
FTC_GEX1 <- RunUMAP(FTC_GEX1, reduction = 'mnn', dims = 1:num_pc, verbose = F)
DimPlot(FTC_GEX1, group.by = 'SampleID')
DimPlot(FTC_GEX1, group.by = 'tissue_type')

#FeaturePlot
FeaturePlot(FTC_GEX1, features = c("TG", "TPO", "TFF3", "VWF", "ACTA2", "LYZ", "CD3D"), combine = T)
FeaturePlot(FTC_GEX1, features = c('CD79A', 'CD19', 'SDC1', 'MZB1', 'CD3E', 'CD3G', 'CD14', 'FCGR3A', 'KIT', 'IL1RL1', 'MS4A2'), ncol = 4)

FeaturePlot(FTC_GEX1, features = c("nCount_RNA", "nFeature_RNA", 'percent.mt'), max.cutoff = 700)

#Reclustering: clustering하기 전 FTC_GEX1으로 다시 clustering
DT::datatable(FTC_GEX1@meta.data)
FTC_GEX1@meta.data$orig.ident
num_PC = 39
FTC_GEX1 <- RunUMAP(FTC_GEX1, reduction = 'mnn', dims = 1:num_PC, verbose = F)

FTC_GEX1 <- FindNeighbors(FTC_GEX1, reduction = 'mnn', dims = 1:num_PC)
FTC_GEX1 <- FindClusters(FTC_GEX1, resolution = 1.0)

saveRDS(FTC_GEX1, '/home/jurny0402/FTC/code/saved_codes/reclustering/230502_FTC_GEX_reclustered_res1.0.rds')
FTC_GEX1 <- readRDS('/home/jurny0402/FTC/code/saved_codes/reclustering/230502_FTC_GEX_reclustered_res1.0.rds')
#중간저장

#Visualizing
DimPlot(FTC_GEX1, label = T, label.size = 4)
DimPlot(FTC_GEX1, group.by = 'SampleID', label = T, label.size = 4, repel = T)
DimPlot(FTC_GEX1, split.by = "tissue_type", group.by = "SampleID", label=T, order = c("Thy1", "Thy4", "Thy5", "Thy6", "Thy7", "Thy10","FT1", "FT2", "FT3", "FT4", "FT5"))
DimPlot(FTC_GEX1, group.by = 'tissue_type', label = F, label.size = 4, repel = T)

table(FTC_GEX1@active.ident)
FTC_GEX1@active.ident <- subset(FTC_GEX1@active.ident, FTC_GEX1@active.ident==0:25) 
FTC_GEX1$seurat_clusters[[57]]

table(FTC_GEX1$seurat_clusters)
FTC_GEX1_filtered <- subset(FTC_GEX1@meta.data, "seurat_clusters" = c(1:26)) # filtered cluster
rownames(FTC_GEX1_filtered[1,"seurat_clusters"])
nrow(FTC_GEX1_filtered)
View(FTC_GEX1$seurat_clusters)

#0:26 cluster만 골라내기 -> seuratobject로 만들기
FTC_GEX1_final <- WhichCells(FTC_GEX1, ident = 0:26)
FTC_GEX2 <- subset(FTC_GEX1, cells=FTC_GEX1_final)

DimPlot(FTC_GEX2, label = T, label.size = 4)

#Finding markers
FTC_GEX2 <- PrepSCTFindMarkers(FTC_GEX2)
saveRDS(FTC_GEX2, '/home/jurny0402/FTC/code/saved_codes/reclustering/230509_FTC_GEX2_prepSCTFindmarkers.rds')
FTC_GEX2 <- readRDS('/home/jurny0402/FTC/code/saved_codes/reclustering/230509_FTC_GEX2_prepSCTFindmarkers.rds')

#clustertree -> 실패
FTC_GEX3 <- CreateSeuratObject(counts = GetAssayData(FTC_GEX2, assay = "RNA", slot = "counts"), meta.data = FTC_GEX2@meta.data)
FTC_GEX3 <- SCTransform(FTC_GEX3)
Idents(FTC_GEX3) <- Idents(FTC_GEX2)
BuildClusterTree(FTC_GEX3, slot = "scale.data")
PlotClusterTree(FTC_GEX3)

saveRDS(FTC_GEX3, file = '/home/jurny0402/FTC/code/saved_codes/reclustering/230509_FTC_GEX3.rds')
FTC_GEX3 <- readRDS('/home/jurny0402/FTC/code/saved_codes/reclustering/230509_FTC_GEX3.rds')

allmarkers <- FindAllMarkers(FTC_GEX2, min.diff.pct = 0.1, only.pos = TRUE, logfc.threshold = 1.0, test.use = 'bimod')
write.csv(allmarkers, file = '/home/jurny0402/FTC/code/saved_codes/reclustering/230509_FTC_GEX2_FTC_GEX2_allmarkers_threshold=1.0.csv')
saveRDS(allmarkers, file = '/home/jurny0402/FTC/code/saved_codes/reclustering/230509_FTC_GEX2_allmarkers_threshold=1.0.rds')
allmarkers <- readRDS('/home/jurny0402/FTC/code/saved_codes/reclustering/230509_FTC_GEX2_allmarkers_threshold=1.0.rds')
#중간저장

allmarkers$cluster
View(allmarkers)

DimPlot(FTC_GEX2, label = T, repel = T)
DimPlot(FTC_GEX2, label = T, repel = T, split.by = 'tissue_type', ncol = 2)
DimPlot(FTC_GEX2, label = T, repel = T, split.by = 'SampleID', ncol = 2)
DimPlot(FTC_GEX2, split.by = "tissue_type", group.by = "SampleID", label=T, order = c("Thy1", "Thy4", "Thy5", "Thy6", "Thy7", "Thy10","FT1", "FT2", "FT3", "FT4", "FT5"))
DimPlot(FTC_GEX2, group.by = 'tissue_type', label = F, label.size = 4, repel = T)

FeaturePlot(FTC_GEX2, features = c('CDH5', 'CLEC14A', 'DIPK2B', 'ERG', 'MMRN2','ROBO4'), combine = T)
FeaturePlot(FTC_GEX2, features = c('ADAMTS3', 'BRINP3', 'ECM2', 'C1QTNF7', 'EPB41L2', 'OLFML1','OMD','SLIT2','SLIT3'), combine = T)
FeaturePlot(FTC_GEX2, features = c('AIF1', 'ALOX5AP', 'C1QA', 'C3AR1', 'CD14', 'FCER1G', 'FOLR2', 'IGSF21', 'RNASE1','TYROBP'), combine = T)
FeaturePlot(FTC_GEX2, features = c('ACTA2', 'ACTG2', 'KCNMA1', 'MYH11', 'PLN', 'PPP1R14A'), combine = T)
FeaturePlot(FTC_GEX2, features = c('CSF3R', 'FCAR', 'FCGR3B', 'RETN', 'S100A12', 'S100A8', 'S100A9'), combine = T)
FeaturePlot(FTC_GEX2, features = c('ACAP1', 'CARD11', 'CCL5', 'CD247', 'CD52', 'CXCR4', 'GPR174','IL2RB','LCK','PRF1','RASGRP1','SAMD3','TRBC1'), combine = T)
FeaturePlot(FTC_GEX2, features = c('ABCC', 'ACP1', 'ADK', 'AHCYL1', 'APLP2', 'APP', 'ARHGAP','IYD','MAOA'), combine = T)
FeaturePlot(FTC_GEX2, features = c('IGJ', 'CD79A', 'DERL3', 'BANK1', 'MS4A1', 'LY9'), combine = T)
FeaturePlot(FTC_GEX2, features=c('KLRK','NCAM1','KLRD','CD69','NCR1','GZMB'))
FeaturePlot(FTC_GEX2, features=c('ITGAM','CD14','CD68','CSF1R','CD163','MRC1'))
FeaturePlot(FTC_GEX2, features = c('IGKV3-23', 'IGHA1', 'IGHG1', 'IGHV1-18','IGHV1-2','IGHV1-3','IGHV1-69D','IGKV1-5','IGKV3-23'), combine = T)

#annotation
FTC_GEX2 <- RenameIdents(FTC_GEX2, '2' = 'Endo', '11' = 'Endo','13' = 'Endo','19' = 'Endo','26' = 'Endo')
FTC_GEX2 <- RenameIdents(FTC_GEX2, '14' = 'Fib','16' = 'sMuscle','25' = 'Neutro')
FTC_GEX2 <- RenameIdents(FTC_GEX2, '3'= 'Macro/TAMs','22'= 'Macro/TAMs','18'= 'NK','23'= 'NK')
FTC_GEX2 <- RenameIdents(FTC_GEX2, '0'= 'T', '6' = 'T','7' = 'T','8' = 'T','9' = 'T','21' = 'T','5' = 'B','12' = 'B','4' = 'Plasma B','24' = 'Plasma B')
FTC_GEX2 <- RenameIdents(FTC_GEX2, '1'= 'Fol','10'= 'Fol','15'= 'Fol','17'= 'Fol','20'= 'Fol')

DimPlot(FTC_GEX2, label = TRUE, label.size = 4)
DimPlot(FTC_GEX2, label = TRUE, label.size = 4, split.by = "tissue_type")

#WhichCells
#clust1 <- WhichCells(FTC_GEX1, idents = "Endo")
#clust2 <- WhichCells(FTC_GEX1, idents = "Fol")
#clust3 <- WhichCells(FTC_GEX1, idents = "Fib")
#clust4 <- WhichCells(FTC_GEX1, idents = "T")
#clust5 <- WhichCells(FTC_GEX1, idents = "B")
#clust6 <- WhichCells(FTC_GEX1, idents = "Macro")
#clust7 <- WhichCells(FTC_GEX1, idents = "Plasma B")
#clust8 <- WhichCells(FTC_GEX1, idents = "sMuscle")
#clust9 <- WhichCells(FTC_GEX1, idents = "Neutro")
#clust10 <- WhichCells(FTC_GEX1, idents = "NK")
#clust11 <- WhichCells(FTC_GEX1, idents = "TAMs")
#저장
saveRDS(clust1, file = '/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=Endo.rds')
saveRDS(clust2, file = '/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=Fol.rds')
saveRDS(clust3, file = '/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=Fib.rds')
saveRDS(clust4, file = '/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=T.rds')
saveRDS(clust5, file = '/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=B.rds')
saveRDS(clust6, file = '/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=Macro.rds')
saveRDS(clust7, file = '/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=Plasma B.rds')
saveRDS(clust8, file = '/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=sMuscle.rds')
saveRDS(clust9, file = '/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=Neutro.rds')
saveRDS(clust10, file = '/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=NK.rds')
saveRDS(clust11, file = '/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=TAMs.rds')

#macro/plasma B 이상
clust1 <- readRDS('/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=Endo.rds')
clust2 <- readRDS('/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=Fol.rds')
clust3 <- readRDS('/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=Fib.rds')
clust4 <- readRDS('/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=T.rds')
clust5 <- readRDS('/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=B.rds')
clust6 <- readRDS('/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=Macro.rds')
clust7 <- readRDS('/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=Plasma B.rds')
clust8 <- readRDS('/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=sMuscle.rds')
clust9 <- readRDS('/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=Neutro.rds')
clust10 <- readRDS('/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=NK.rds')
clust11 <- readRDS('/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=TAMs.rds')

DimPlot(FTC_GEX1, label = TRUE, label.size = 4, cells.highlight = clust1, sizes.highlight = .1, split.by = "tissue_type")

---------------------------------------------------------------------------------
z <- subset(FTC_GEX1@meta.data, "SampleID" = "FT1")
nrow(FTC_GEX1@meta.data)
View(FTC_GEX1@meta.data)
FTC_GEX1$SampleID[[2]]
FTC_GEX1$seurat_clusters[[10000]]

FTC_GEX1@meta.data[[2,"seurat_clusters"]]

FT1=vector();FT2=vector();FT3=vector();FT4=vector();FT5=vector()
Thy1=vector();Thy10=vector();Thy4=vector();Thy5=vector();Thy6=vector();Thy7=vector();
j=1
for (i in 1:73619) {
    if (FTC_GEX1@meta.data[[i,"SampleID"]]=="FT1") {
      FT1[j] <- c(FTC_GEX1@meta.data[[i,"seurat_clusters"]])
      j <- j+1
  }
}
sort(FT1)
length(FT1)
sort(unique(FT1))

j=1
for (i in 1:73619) {
  if (FTC_GEX1@meta.data[[i,"SampleID"]]=="FT2") {
    FT2[j] <- c(FTC_GEX1@meta.data[[i,"seurat_clusters"]])
    j <- j+1
  }
}
sort(unique(FT2))

j=1
for (i in 1:73619) {
  if (FTC_GEX1@meta.data[[i,"SampleID"]]=="FT3") {
    FT3[j] <- c(FTC_GEX1@meta.data[[i,"seurat_clusters"]])
    j <- j+1
  }
}
sort(unique(FT3))

j=1
for (i in 1:73619) {
  if (FTC_GEX1@meta.data[[i,"SampleID"]]=="FT4") {
    FT4[j] <-  c(FTC_GEX1@meta.data[[i,"seurat_clusters"]])
    j <- j+1
  }
}
sort(unique(FT4))

j=1
for (i in 1:73619) {
  if (FTC_GEX1@meta.data[[i,"SampleID"]]=="FT5") {
    FT5[j] <-  c(FTC_GEX1@meta.data[[i,"seurat_clusters"]])
    j <- j+1
  }
}
sort(unique(FT5))

j=1
for (i in 1:73619) {
  if (FTC_GEX1@meta.data[[i,"SampleID"]]=="Thy1") {
    Thy1[j] <-  c(FTC_GEX1@meta.data[[i,"seurat_clusters"]])
    j <- j+1
  }
}
sort(Thy1)
levels(Thy1)

j=1
for (i in 1:73619) {
  if (FTC_GEX1@meta.data[[i,"SampleID"]]=="FT1") {
    FT1[j] <-  c(FTC_GEX1@meta.data[[i,"seurat_clusters"]])
    j <- j+1
  }
}
sort(FT1)
levels(FT1)

j=1
for (i in 1:73619) {
  if (FTC_GEX1@meta.data[[i,"SampleID"]]=="Thy10") {
    Thy10[j] <-  c(FTC_GEX1@meta.data[[i,"seurat_clusters"]])
    j <- j+1
  }
}
sort(Thy10)
levels(Thy10)

j=1
for (i in 1:73619) {
  if (FTC_GEX1@meta.data[[i,"SampleID"]]=="Thy4") {
    Thy4[j] <-  c(FTC_GEX1@meta.data[[i,"seurat_clusters"]])
    j <- j+1
  }
}
sort(Thy4)
levels(Thy4)

j=1
for (i in 1:73619) {
  if (FTC_GEX1@meta.data[[i,"SampleID"]]=="Thy5") {
    Thy5[j] <-  c(FTC_GEX1@meta.data[[i,"seurat_clusters"]])
    j <- j+1
  }
}
sort(Thy5)
levels(Thy5)

j=1
for (i in 1:73619) {
  if (FTC_GEX1@meta.data[[i,"SampleID"]]=="Thy6") {
    Thy6[j] <-  c(FTC_GEX1@meta.data[[i,"seurat_clusters"]])
    j <- j+1
  }
}
sort(Thy6)
levels(Thy6)

j=1
for (i in 1:73619) {
  if (FTC_GEX1@meta.data[[i,"SampleID"]]=="Thy7") {
    Thy7[j] <-  c(FTC_GEX1@meta.data[[i,"seurat_clusters"]])
    j <- j+1
  }
}
sort(Thy7)
levels(Thy7)
------------------------------------------------------------------
FT1 <- c(0:31,33:40,42,44:45,47:49,53);FT2 <- c(0:3,5:28,30:43,45:46,48:51,53,55);FT3 <- c(0:46,48:57);FT4 <- c(0:3, 5:23,25,27,30:31,35,36,39:42,44,46,50,55);FT5 <- c(0:3,5:23,25,27,30:33,35:36,38:40,42);Thy1 <- c(0:3,5:27,29:30,35:36,38:40,42,44,46);Thy4 <- c(0:3,5:16,18:23,25:27,30,32,36,39:42);Thy5 <- c(0:3,5:16,18:27,29:31,35:36,38,40,42);Thy6 <- c(0:1,3,5:16,18:23,25:27,29:30,36,38,40);Thy7 <- c(1,15,17,22,41,46);Thy10 <- c(0:3,5:16,18:27,29:30,36,39:40,42,44,46)

FT1_exists=vector();FT2_exists=vector();FT3_exists=vector();FT4_exists=vector();FT5_exists=vector()
Thy1_exists=vector();Thy4_exists=vector();Thy5_exists=vector();Thy6_exists=vector();Thy7_exists=vector();Thy10_exists=vector()
number <- c(0:57)
FT1_exists <- c(0:57)
for (i in 1:length(FT1)) {
  for (j in 1:58) {
    if (FT1[i] == FT1_exists[j]) {
      FT1_exists[j] <- c("o")
    }
  }
}

FT2_exists <- c(0:57)
for (i in 1:length(FT2)) {
  for (j in 1:58) {
    if (FT2[i] == FT2_exists[j]) {
      FT2_exists[j] <- c("o")
    }
  }
}

FT3_exists <- c(0:57)
for (i in 1:length(FT3)) {
  for (j in 1:58) {
    if (FT3[i] == FT3_exists[j]) {
      FT3_exists[j] <- c("o")
    }
  }
}

FT4_exists <- c(0:57)
for (i in 1:length(FT4)) {
  for (j in 1:58) {
    if (FT4[i] == FT4_exists[j]) {
      FT4_exists[j] <- c("o")
    }
  }
}

FT5_exists <- c(0:57)
for (i in 1:length(FT5)) {
  for (j in 1:58) {
    if (FT5[i] == FT5_exists[j]) {
      FT5_exists[j] <- c("o")
    }
  }
}

Thy1_exists <- c(0:57)
for (i in 1:length(Thy1)) {
  for (j in 1:58) {
    if (Thy1[i] == Thy1_exists[j]) {
      Thy1_exists[j] <- c("o")
    }
  }
}

Thy4_exists <- c(0:57)
for (i in 1:length(Thy4)) {
  for (j in 1:58) {
    if (Thy4[i] == Thy4_exists[j]) {
      Thy4_exists[j] <- c("o")
    }
  }
}

Thy5_exists <- c(0:57)
for (i in 1:length(Thy5)) {
  for (j in 1:58) {
    if (Thy5[i] == Thy5_exists[j]) {
      Thy5_exists[j] <- c("o")
    }
  }
}

Thy6_exists <- c(0:57)
for (i in 1:length(Thy6)) {
  for (j in 1:58) {
    if (Thy6[i] == Thy6_exists[j]) {
      Thy6_exists[j] <- c("o")
    }
  }
}

Thy7_exists <- c(0:57)
for (i in 1:length(Thy7)) {
  for (j in 1:58) {
    if (Thy7[i] == Thy7_exists[j]) {
      Thy7_exists[j] <- c("o")
    }
  }
}

Thy10_exists<- c(0:57)
for (i in 1:length(Thy10)) {
  for (j in 1:58) {
    if (Thy10[i] == Thy10_exists[j]) {
      Thy10_exists[j] <- c("o")
    }
  }
}

FT <- cbind(number,FT1_exists)
FT <- cbind(FT,FT2_exists)
FT <- cbind(FT,FT3_exists)
FT <- cbind(FT,FT4_exists)
FT <- cbind(FT,FT5_exists)
Thy <- cbind(number,Thy1_exists)
Thy <- cbind(Thy, Thy4_exists)
Thy <- cbind(Thy, Thy5_exists)
Thy <- cbind(Thy, Thy6_exists)
Thy <- cbind(Thy, Thy7_exists)
Thy <- cbind(Thy, Thy10_exists)
cluster <- cbind(FT,Thy)
View(cluster)

write.csv(cluster, file = '/home/jurny0402/FTC/code/saved_codes/reclustering/230523cluster.csv')

#cluster 골라내기-Option 1&2
FTC_GEX1_final <- WhichCells(FTC_GEX1, ident = c(0:3,5:16,18:27,29:30,36,38:40,42))
FTC_GEX2 <- subset(FTC_GEX1, cells=FTC_GEX1_final)

DimPlot(FTC_GEX2, label = T, label.size = 4)
DimPlot(FTC_GEX2, group.by = 'tissue_type')

#annotation
FeaturePlot(FTC_GEX2, features = c('CDH5', 'CLEC14A', 'DIPK2B', 'ERG', 'MMRN2','ROBO4'), combine = T)
FeaturePlot(FTC_GEX2, features = c('ADAMTS3', 'BRINP3', 'ECM2', 'C1QTNF7', 'EPB41L2', 'OLFML1','OMD','SLIT2','SLIT3'), combine = T)
FeaturePlot(FTC_GEX2, features = c('AIF1', 'ALOX5AP', 'C1QA', 'C3AR1', 'CD14', 'FCER1G', 'FOLR2', 'IGSF21', 'RNASE1','TYROBP'), combine = T)
FeaturePlot(FTC_GEX2, features = c('ACTA2', 'ACTG2', 'KCNMA1', 'MYH11', 'PLN', 'PPP1R14A'), combine = T)
FeaturePlot(FTC_GEX2, features = c('CSF3R', 'FCAR', 'FCGR3B', 'RETN', 'S100A12', 'S100A8', 'S100A9'), combine = T)
FeaturePlot(FTC_GEX2, features = c('ACAP1', 'CARD11', 'CCL5', 'CD247', 'CD52', 'CXCR4', 'GPR174','IL2RB','LCK','PRF1','RASGRP1','SAMD3','TRBC1'), combine = T)
FeaturePlot(FTC_GEX2, features = c('ABCC', 'ACP1', 'ADK', 'AHCYL1', 'APLP2', 'APP', 'ARHGAP','IYD','MAOA'), combine = T)
FeaturePlot(FTC_GEX2, features = c('IGJ', 'CD79A', 'DERL3', 'BANK1', 'MS4A1', 'LY9'), combine = T)
FeaturePlot(FTC_GEX2, features=c('KLRK','NCAM1','KLRD','CD69','NCR1','GZMB'))
FeaturePlot(FTC_GEX2, features=c('ITGAM','CD14','CD68','CSF1R','CD163','MRC1'))
FeaturePlot(FTC_GEX2, features = c('IGKV3-23', 'IGHA1', 'IGHG1', 'IGHV1-18','IGHV1-2','IGHV1-3','IGHV1-69D','IGKV1-5','IGKV3-23'), combine = T)

#annotation
FTC_GEX2 <- RenameIdents(FTC_GEX2, '2' = 'Endo', '11' = 'Endo','13' = 'Endo','19' = 'Endo','26' = 'Endo','27' = 'Endo','39' = 'Endo','40' = 'Endo')
FTC_GEX2 <- RenameIdents(FTC_GEX2, '14' = 'Fib','16' = 'sMuscle','36' = 'sMuscle','42' = 'sMuscle','25' = 'Neutro')
FTC_GEX2 <- RenameIdents(FTC_GEX2, '3'= 'Macro/TAMs','22'= 'Macro/TAMs','18'= 'NK','23'= 'NK')
FTC_GEX2 <- RenameIdents(FTC_GEX2, '0'= 'T', '6' = 'T','7' = 'T','8' = 'T','9' = 'T','21' = 'T','5' = 'B','12' = 'B','29' = 'Plasma B','24' = 'Plasma B')
FTC_GEX2 <- RenameIdents(FTC_GEX2, '1'= 'Fol','10'= 'Fol','15'= 'Fol','20'= 'Fol','30'='TAMs','38'='TAMs')
DimPlot(FTC_GEX2, label = TRUE, label.size = 4)
DimPlot(FTC_GEX2, label = TRUE, label.size = 4, split.by = "tissue_type")

clust1 <- readRDS('/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=Endo.rds')
clust2 <- readRDS('/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=Fol.rds')
clust3 <- readRDS('/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=Fib.rds')
clust4 <- readRDS('/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=T.rds')
clust5 <- readRDS('/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=B.rds')
clust6 <- readRDS('/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=Macro.rds')
clust7 <- readRDS('/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=Plasma B.rds')
clust8 <- readRDS('/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=sMuscle.rds')
clust9 <- readRDS('/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=Neutro.rds')
clust10 <- readRDS('/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=NK.rds')
clust11 <- readRDS('/home/jurny0402/FTC/code/saved_codes/WhichCells/230502_FTC_GEX1_WhichCells=TAMs.rds')

DimPlot(FTC_GEX2, label = TRUE, label.size = 4, cells.highlight = clust11, sizes.highlight = .1, split.by = "tissue_type")
---------------------------------------------------------------------------------------------------------------------------------
#Option 3
#cluster 골라내기
FTC_GEX1_final <- WhichCells(FTC_GEX1, ident = c(0:3,5:16,18:23,25:27,29:30,36,40,42))
FTC_GEX2 <- subset(FTC_GEX1, cells=FTC_GEX1_final)

DimPlot(FTC_GEX2, label = T, label.size = 4)
DimPlot(FTC_GEX2, group.by = 'tissue_type')
  
#annotation
FTC_GEX2 <- RenameIdents(FTC_GEX2, '2' = 'Endo', '11' = 'Endo','13' = 'Endo','19' = 'Endo','26' = 'Endo','27' = 'Endo','40' = 'Endo')
FTC_GEX2 <- RenameIdents(FTC_GEX2, '14' = 'Fib','16' = 'sMuscle','36' = 'sMuscle','42' = 'sMuscle','25' = 'Neutro')
FTC_GEX2 <- RenameIdents(FTC_GEX2, '3'= 'Macro/TAMs','22'= 'Macro/TAMs','18'= 'NK','23'= 'NK')
FTC_GEX2 <- RenameIdents(FTC_GEX2, '0'= 'T', '6' = 'T','7' = 'T','8' = 'T','9' = 'T','21' = 'T','5' = 'B','12' = 'B','29' = 'Plasma B')
FTC_GEX2 <- RenameIdents(FTC_GEX2, '1'= 'Fol','10'= 'Fol','15'= 'Fol','20'= 'Fol','30'='TAMs')
DimPlot(FTC_GEX2, label = TRUE, label.size = 4)
DimPlot(FTC_GEX2, label = TRUE, label.size = 4, split.by = "tissue_type")
  
FTC_GEX2_final <- WhichCells(FTC_GEX2, ident = "T") 
FTC_GEX3 <- subset(FTC_GEX2, cells = FTC_GEX2_final)
DimPlot(FTC_GEX3, label = TRUE, label.size = 4)
saveRDS(FTC_GEX3, '/home/jurny0402/FTC/code/saved_codes/reclustering/230704_FTC_GEX3_onlyTclusters.rds')
FTC_GEX3 <- readRDS('/home/jurny0402/FTC/code/saved_codes/reclustering/230704_FTC_GEX3_onlyTclusters.rds')

FTC_GEX2_final <- WhichCells(FTC_GEX2, ident = "Endo") 
FTC_GEX3 <- subset(FTC_GEX2, cells = FTC_GEX2_final)

DotPlot(FTC_GEX3, features = c("VEGFA","PDGFA","PDGFB","VEGFR1","VEGFR2","eNOS"," VEGFR3"),split.by = "tissue_type",scale = TRUE)
DotPlot(FTC_GEX2, features = c("VEGFA","PDGFA","PDGFB","VEGFR1","VEGFR2","eNOS"," VEGFR3"),scale = TRUE)

View(FTC_GEX3@meta.data$)
RidgePlot(FTC_GEX2, features = c("VEGFA","PDGFA","PDGFB","VEGFR1","VEGFR2","eNOS"," VEGFR3"), ncol = 2)
FeaturePlot(FTC_GEX2, features = c("VEGFA","PDGFA","PDGFB","VEGFR1","VEGFR2","eNOS"," VEGFR3"),split.by = "tissue_type")
