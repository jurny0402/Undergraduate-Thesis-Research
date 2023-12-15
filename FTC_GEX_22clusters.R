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

#Clustering
num_PC = 39
FTC_GEX1 <- RunUMAP(FTC_GEX1, reduction = 'mnn', dims = 1:39, verbose = F)

FTC_GEX1 <- FindNeighbors(FTC_GEX1, reduction = 'mnn', dims = 1:13)

FTC_GEX1 <- FindClusters(FTC_GEX1, resolution = 0.3)

saveRDS(FTC_GEX1, '/home/jurny0402/FTC/code/saved_codes/230408_FTC_GEX_clustered_res0.3_dims1:13.rds')
FTC_GEX1 <- readRDS('/home/jurny0402/FTC/code/saved_codes/230408_FTC_GEX_clustered_res0.3_dims1:13.rds')

#중간저장

#Visualizing
DimPlot(FTC_GEX1, label = T, label.size = 4)
DimPlot(FTC_GEX1, group.by = 'SampleID', label = T, label.size = 4, repel = T)
DimPlot(FTC_GEX1, split.by = "tissue_type", group.by = "SampleID", label=T, order = c("Thy1", "Thy4", "Thy5", "Thy6", "Thy7", "Thy10","FT1", "FT2", "FT3", "FT4", "FT5"))

##clustering, visualizing 과정을 dims = 1:40로 다시 해보기
table(FTC_GEX1@active.ident)

#Finding markers
FTC_GEX1 <- PrepSCTFindMarkers(FTC_GEX1)
saveRDS(FTC_GEX1, '/home/jurny0402/FTC/code/saved_codes/230425_FTC_GEX_22cl_prepSCTFindmarkers.rds')

FTC_GEX1 <- readRDS('/home/jurny0402/FTC/code/saved_codes/230425_FTC_GEX_22cl_prepSCTFindmarkers.rds')
#중간저장

allmarkers <- FindAllMarkers(FTC_GEX1, min.diff.pct = 0.1, only.pos = TRUE, logfc.threshold = 0.25, test.use = 'bimod')
write.csv(allmarkers, file = '/home/jurny0402/FTC/code/saved_codes/23042_FTC_GEX1_22cl_allmarkers.csv')
saveRDS(allmarkers, file = '/home/jurny0402/FTC/code/saved_codes/230425_FTC_GEX1_22cl_allmarkers.rds')
allmarkers <- readRDS('/home/jurny0402/FTC/code/saved_codes/230425_FTC_GEX1_22cl_allmarkers.rds')
#중간저장

allmarkers$cluster

DimPlot(FTC_GEX1, label = T, repel = T)
DimPlot(FTC_GEX1, label = T, repel = T, split.by = 'tissue_type', ncol = 1)
DimPlot(FTC_GEX1, label = T, repel = T, split.by = 'SampleID', ncol = 1)

FeaturePlot(FTC_GEX1, features = c('TPO', 'TG', 'KRT18', 'KRT19', 'TSHR'), combine = T)
FeaturePlot(FTC_GEX1, features = c('RAMP2', 'VWF', 'ARL15', 'PLPP1', 'PTPRG', 'STC1'), combine = T)
FeaturePlot(FTC_GEX1, features = c('PDGFRA', 'DCN', 'COL3A1', 'COL1A1', 'MGP', 'COL1A2', 'ACTA2', 'FAP', 'TAGLN'), combine = T)
FeaturePlot(FTC_GEX1, features = c('CD3D', 'CD3E', 'FOXP3'), combine = T)
FeaturePlot(FTC_GEX1, features = c('IGJ', 'CD79A', 'DERL3', 'BANK1', 'MS4A1', 'LY9'), combine = T)
FeaturePlot(FTC_GEX1, features = c('CD68', 'IL1B', 'FCGR3A', 'TNF', 'S100A9', 'CD14', 'CTSD', 'CD163'), combine = T)

#jurny's annotation
FeaturePlot(FTC_GEX1, features = c('CDH5', 'CLEC14A', 'DIPK2B', 'ERG', 'MMRN2','ROBO4'), combine = T)
FeaturePlot(FTC_GEX1, features = c('ADAMTS3', 'BRINP3', 'ECM2', 'C1QTNF7', 'EPB41L2', 'OLFML1','OMD','SLIT2','SLIT3'), combine = T)
FeaturePlot(FTC_GEX1, features = c('AIF1', 'ALOX5AP', 'C1QA', 'C3AR1', 'CD14', 'FCER1G', 'FOLR2', 'IGSF21', 'RNASE1','TYROBP'), combine = T)
FeaturePlot(FTC_GEX1, features = c('IGHA1', 'IGHG1', 'IGHG3','IGHV1-18','IGHV1-3','IGHV3-30','IGHV4-39','IGKC'), combine = T)
FeaturePlot(FTC_GEX1, features = c('ACTA2', 'ACTG2', 'KCNMA1', 'MYH11', 'PLN', 'PPP1R14A'), combine = T)
FeaturePlot(FTC_GEX1, features = c('CSF3R', 'FCAR', 'FCGR3B', 'RETN', 'S100A12', 'S100A8', 'S100A9'), combine = T)
FeaturePlot(FTC_GEX1, features = c('ACAP1', 'CARD11', 'CCL5', 'CD247', 'CD52', 'CXCR4', 'GPR174','IL2RB','LCK','PRF1','RASGRP1','SAMD3','TRBC1'), combine = T)
FeaturePlot(FTC_GEX1, features = c('ABCC', 'ACP1', 'ADK', 'AHCYL1', 'APLP2', 'APP', 'ARHGAP','IYD','MAOA'), combine = T)
FeaturePlot(FTC_GEX1, features = c('IGJ', 'CD79A', 'DERL3', 'BANK1', 'MS4A1', 'LY9'), combine = T)
FeaturePlot(FTC_GEX1, features=c('GZMA','NKG7'))
FeaturePlot(FTC_GEX1, features=c('KLRK','NCAM1','KLRD','CD69','NCR1','GZMB'))
FeaturePlot(FTC_GEX1, features=c('ITGAM','CD14','CD68','CSF1R','CD163','MRC1'))
FeaturePlot(FTC_GEX1, features=c('CD3D','CD8A','CD4','CD19','NCAM1','PTPRC'))

#Cell type annotation
FTC_GEX1 <- RenameIdents(FTC_GEX1, '3' = 'Fol', '6' = 'Fol', '12' = 'Mac','15' = 'Fol','16' = 'Fol')
FTC_GEX1 <- RenameIdents(FTC_GEX1, '2' = 'Endo', '10' = 'Endo','13'='T/ILC','14'='Endo','17'='Endo')
FTC_GEX1 <- RenameIdents(FTC_GEX1, '8' = 'Fib', '11' = 'Fib','21'='Fib')
FTC_GEX1 <- RenameIdents(FTC_GEX1, '0' = 'T/ILC', '1' = 'NK')
FTC_GEX1 <- RenameIdents(FTC_GEX1, '4' = 'B','20'='B','18'='B')
FTC_GEX1 <- RenameIdents(FTC_GEX1, '5' = 'TAMs','19'='Mac','9'='Mac','7'='Mac')

#jurny's cell type annotation
FTC_GEX1 <- RenameIdents(FTC_GEX1, '2' = 'Endo', '10' = 'Endo','14'='Endo','17'='Endo')
FTC_GEX1 <- RenameIdents(FTC_GEX1, '11' = 'Fib','8'="sMuscle")
FTC_GEX1 <- RenameIdents(FTC_GEX1, '5'= 'TAMs','9'='Plasma B','12'='Macro','7'= 'Plasma B','18'='Plasma B','21'='Plasma B')
FTC_GEX1 <- RenameIdents(FTC_GEX1, '0' = 'T', '1'= 'NK', '13' = 'T','19'='Neutro')
FTC_GEX1 <- RenameIdents(FTC_GEX1, '3'= 'Fol','6'= 'Fol','15'= 'Fol','16'='Fol')
FTC_GEX1 <- RenameIdents(FTC_GEX1, '4' = 'B','20'='B')

DimPlot(FTC_GEX1, label = TRUE, label.size = 4)
DimPlot(FTC_GEX1, label = TRUE, label.size = 4, split.by = "tissue_type")

#workspace에서 삭제
rm(a);rm(a1);rm(a2);rm(b);rm(b1);rm(b2);rm(c);rm(c1);rm(c2);rm(clust1);rm(d);rm(d1);rm(d2)
rm(f);rm(f1);rm(f2);rm(f3);rm(f4);rm(t);rm(m_a);rm(m_b);rm(m_c);rm(m_d);rm(m_f)
rm(FTC_GEX3);rm(FTC_GEX4);rm(FTC_GEX2);rm(FTC_GEX_BR)
rm(meta_mat);rm(meta.info);rm(degs_list1)
rm(batch_list);rm(ident,ident2,j,n,o)
#garbage reset
gc(reset=T)

view(allmarkers)

#macrophages
a <- c(6,8,14,15,16,21,25,28,29,32,6,8,14,15,21,23,24,25,26,28,29,30,31,33,6,8,12,16,21,25,33,
       6,8,15,25,26,28,29,31,33,21,6,14,15,16,21,25,28,32,33,6,13,16,19,25,28,29,32,6,16,29)

a1 <- sort.int(a);n=0
for (i in 1:length(a1)) {
  if (a1[i] == a1[i+1]) {
    a1[i]=0
  }
  if (a1[i] != 0) {
    n <- n+1 
  }
}

a1 <- a1[!a1 %in% 0]
length(a1)
levels(a1)

a2 <- rep(0, times=length(a1))

for (i in 1:length(a)) {
  for (j in 1:length(a1)) {
    if (a[i] == a1[j]) {
      a2[j] <- a2[j]+1
    }
  }
}

m_a <- rbind(a1,a2)
m_a <- cbind(a1,a2)
view(m_a)

b <- c(2,10,17,11,5,7,9,12,19,7,9,18,21,8,19,0,1,13,3,6,15,16,4)
b1 <- sort.int(b);n=0
for (i in 1:length(b1)) {
  if (b1[i] == b1[i+1]) {
    b1[i]=0
  }
  if (b1[i] != 0) {
    n <- n+1 
  }
}

b1 <- b1[!b1 %in% 0]
length(b1)
levels(b1)

b2 <- rep(0, times=length(b1))

for (i in 1:length(b)) {
  for (j in 1:length(b1)) {
    if (b[i] == b1[j]) {
      b2[j] <- b2[j]+1
    }
  }
}

m_b <- rbind(b1,b2)
m_b <- cbind(b1,b2)
view(m_b)
