library(ggplot2)
library(hrbrthemes)
load("/home/jurny0402/231116_nor,FTC_seurat/data/231214_adenoma_pairing.rda")
load("/home/jurny0402/231116_nor,FTC_seurat/data/231214_carcinoma_pairing.rda")
load("/home/jurny0402/231116_nor,FTC_seurat/data/231214_normal_pairing.rda")
ggplot(adenoma, aes(Var1, Var2, fill= 0.05)) + 
  geom_tile() +
  scale_fill_gradient(low="black", high="red")
adenoma
