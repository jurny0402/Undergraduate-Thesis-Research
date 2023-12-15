a = data
a2 = view@meta.data
head(a)
head(a2)
a3 = merge(a,a2,"contig_id")
dim(a)
dim(a2)
dim(a3)

saveRDS(a3, file = '/home/jurny0402/TCR_explore/code/merged_data.rds')
merged_data <- read.csv('/home/jurny0402/TCR_explore/code/merged_data.rds')
