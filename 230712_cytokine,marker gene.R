library(pheatmap)
library(Seurat)
library(DESeq2)

FT_data <- readRDS('/home/jurny0402/TCR_explore/code/T cluster only/230704_onlyT_FTC_GEX4_after_clutering.rds')
DimPlot(FT_data, label=T)
DoHeatmap(FT_data, features = c("TCR7","ID3","CCR7","AQP3","SELL"))
#CD4 rest
FeaturePlot(FT_data, features = c("TCF7","ID3","CCR7","AQP3","SELL"))
DoHeatmap(FT_data, features = c("TCF7","ID3","CCR7","AQP3","SELL"))
#CD4 act1/2
FeaturePlot(FT_data, features = c("IL4R","STAT1","MAL","SOCS1"))
DoHeatmap(FT_data, features = c("IL4R","STAT1","MAL","SOCS1"))
#CD4 act3
FeaturePlot(FT_data, features = c("IL2","ODC1","WARS","PYCR1","TNF","MIR155HG","NME1"))
DoHeatmap(FT_data, features = c("IL2","ODC1","WARS","PYCR1","TNF","MIR155HG","NME1"))
#CD4 Treg
FeaturePlot(FT_data, features = c("FOXP3","RGS1","STAT1","LGALS3","IL2RA","CTLA4","TIGIT","TNFRSF4","TNFRSF18"))
DoHeatmap(FT_data, features = c("FOXP3","RGS1","STAT1","LGALS3","IL2RA","CTLA4","TIGIT","TNFRSF4","TNFRSF18"))
#CD8 EM/TRM rest
FeaturePlot(FT_data, features = c("KLRB1","AMICA1","LGALS1","GZMK","CXCR6","ANXA2","ITGA1"))
DoHeatmap(FT_data, features = c("KLRB1","AMICA1","LGALS1","GZMK","CXCR6","ANXA2","ITGA1"))
#CD8 TRM rest
FeaturePlot(FT_data, features = c("AMICA1","PDZD8","CXCR6","RCAN2","PRF1","ANXA2","CRIP1","ITGA1","GZMB"))
DoHeatmap(FT_data, features = c("AMICA1","PDZD8","CXCR6","RCAN2","PRF1","ANXA2","CRIP1","ITGA1","GZMB"))
#CD8 EM/TRM act
FeaturePlot(FT_data, features = c("GZMB","IFNG","HOPX","CCL4","CCL3","XCL1"))
DoHeatmap(FT_data, features = c("GZMB","IFNG","HOPX","CCL4","CCL3","XCL1"))
#CD8 TEMRA
FeaturePlot(FT_data, features = c("NKG7","GZMK","GNLY","IFIT3","KLRD1","MYO1F","CCL5"))
DoHeatmap(FT_data, features = c("NKG7","GZMK","GNLY","IFIT3","KLRD1","MYO1F","CCL5"))
#CD8 TEMRA rest
FeaturePlot(FT_data, features = c("GZMH","NKG7","MYO1F","CCL5","GNLY"))
DoHeatmap(FT_data, features = c("GZMH","NKG7","MYO1F","CCL5","GNLY"))
#CD8 TEMRA act
FeaturePlot(FT_data, features = c("CCL5","GNLY","IFIT3","KLRD1","XCL2","HOPX","RCAN2","PRF1"))
DoHeatmap(FT_data, features = c("CCL5","GNLY","IFIT3","KLRD1","XCL2","HOPX","RCAN2","PRF1"))

DoHeatmap(FT_data, features = c("TCF7","ID3","CCR7","AQP3","SELL","IL4R","STAT1","MAL","SOCS1","IL2","ODC1","WARS","PYCR1","TNF","MIR155HG","NME1",
                                "FOXP3","RGS1","STAT1","LGALS3","IL2RA","CTLA4","TIGIT","TNFRSF4","TNFRSF18",
                                "KLRB1","AMICA1","LGALS1","GZMK","CXCR6","ANXA2","ITGA1",
                                "AMICA1","PDZD8","CXCR6","RCAN2","PRF1","ANXA2","CRIP1","ITGA1","GZMB",
                                "GZMB","IFNG","HOPX","CCL4","CCL3","XCL1",
                                "NKG7","GZMK","GNLY","IFIT3","KLRD1","MYO1F","CCL5",
                                "GZMH","NKG7","MYO1F","CCL5","GNLY",
                                "CCL5","GNLY","IFIT3","KLRD1","XCL2","HOPX","RCAN2","PRF1"))

FTC <- readRDS('/home/jurny0402/TCR_explore/code/T cluster only/230704_FTC_GEX4.rds')

DoHeatmap(FTC, features = c("TCF7","ID3","CCR7","AQP3","SELL",
                            "IL4R","STAT1","MAL","SOCS1",
                            "IL2","ODC1","WARS","PYCR1","TNF","MIR155HG","NME1",
                                "FOXP3","RGS1","STAT1","LGALS3","IL2RA","CTLA4","TIGIT","TNFRSF4","TNFRSF18",
                                "KLRB1","AMICA1","LGALS1","GZMK","CXCR6","ANXA2","ITGA1",
                                "AMICA1","PDZD8","CXCR6","RCAN2","PRF1","ANXA2","CRIP1","ITGA1","GZMB",
                                "GZMB","IFNG","HOPX","CCL4","CCL3","XCL1",
                                "NKG7","GZMK","GNLY","IFIT3","KLRD1","MYO1F","CCL5",
                                "GZMH","NKG7","MYO1F","CCL5","GNLY",
                                "CCL5","GNLY","IFIT3","KLRD1","XCL2","HOPX","RCAN2","PRF1"))

#Transcriptome profiling of human thymic CD4+ and CD8+ T cells compared to primary peripheral T cells
DoHeatmap(FT_data, features = c("S1PR5","PLEKHG3","TBX21")) #CD8
DoHeatmap(FT_data, features = c("IL6R","IL4R")) #CD4
DoHeatmap(FT_data, features = c("TBX21","IL6R","IL4R"))   

DoHeatmap(FT_data, features = c("CACBA1I","IL6R","IL4R","AK5")) #CD4
DoHeatmap(FT_data, features = c("NABP1","KIF5C","NEO1","RORC","AGAP1","PPP2R2B","PTCH1","SCART1","MIAT")) #CD8
DoHeatmap(FT_data, features = c("CACBA1I","IL6R","IL4R","AK5","NABP1","KIF5C","NEO1","RORC","AGAP1","PPP2R2B","PTCH1","SCART1","MIAT"))

DoHeatmap(FT_data, features = c("ITGB1","CD4","IL6R","IL2RA","IL7R","ITK","ITGA6")) #CD4
DoHeatmap(FT_data, features = c("CCR2","IL10RA","IL32","SELPLG","ITGB7","IL18R1","CCR5","CXCR6","CCL5","NFATC2","IL12RB1","CARD11","ADAM8","RUNX3","CR1",
                                "CXCR5","CCR1","CXCR3","PECAM1")) #CD8
data <- DoHeatmap(FT_data, features = c("CCR2","IL10RA","IL32","SELPLG","ITGB7","IL18R1","CCR5","CXCR6","CCL5","NFATC2","IL12RB1","CARD11","ADAM8","RUNX3","CR1",
                                        "CXCR5","CCR1","CXCR3","PECAM1")) #CD8
View(data$data)
View(FT_data$seurat_clusters)
