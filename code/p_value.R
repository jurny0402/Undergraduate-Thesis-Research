library(ggpubr)
library(rstatix)

df <- readRDS('/home/jurny0402/231116_nor,FTC_seurat/data/231212_cdr3length_distribution.rds')

View(df)
class(df)
colnames(df)
df$values
df$ID

#box plots with p-values
stat.test <- df %>%
  t_test(length ~ ID) %>%
  add_significance()
stat.test

bxp <- ggboxplot(df, x = "ID", y = "length", fill = "#00AFBB")
stat.test <- stat.test %>% add_xy_position(x = "ID")
bxp + 
  stat_pvalue_manual(stat.test, label = "p") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#Grouped data
stat.test <- df %>%
  group_by(dose) %>%
  t_test(length ~ ID) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test

stat.test <- stat.test %>% add_xy_position(x = "values")
bxp <- ggboxplot(df, x = "values", y = "length", fill = "#00AFBB",
                 facet.by = "ID")
bxp + 
  stat_pvalue_manual(stat.test, label = "p.adj") +
  scale_y_continuous(expand = expansion(c("a","c","n")))

#show p-values if significant otherwise show ns
stat.test <- stat.test %>% add_xy_position(x = "ID")
stat.test$custom.label <- ifelse(stat.test$p.adj <= 0.5, stat.test$p.adj, "ns")

bxp + 
  stat_pvalue_manual(stat.test, label = "custom.label") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10)))

#format p-values using the accuracy option
stat.test <- stat.test %>% add_xy_position(x = "ID")
stat.test$p.format <- p_format(
  stat.test$p.adj, accuracy = 0.01,
  leading.zero = FALSE)

bxp + 
  stat_pvalue_manual(stat.test, label = "p.format") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10)))

View(df)

#create plots with significance levels
stat.test <- stat.test %>% add_xy_position(x = "ID")
bxp + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01)
#bar plot
stat.test <- stat.test %>% add_xy_position(fun = "mean_sd", x = "ID")
bxp + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01)

#compare paired samples
stat.test <- df %>%
  t_test(length ~ ID, paired = TRUE) %>%
  add_significance()
stat.test

bxp <-  ggpaired(df, x = "ID", y = "length", fill = "#E7B800",
                 line.color = "gray", line.size = 0.4)
stat.test <- stat.test %>% add_xy_position(x = "ID")
bxp + stat_pvalue_manual(stat.test, label = "p.adj.signif")


stat.test <- stat.test %>% add_xy_position(x = "ID")
bxp + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01)


p <- ggboxplot(df, x = "ID", y = "length",
               color = "ID", palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()
# Change method
p + stat_compare_means(method = "t.test")


my_comparisons <- list( c("a", "c"), c("a", "n"), c("c", "n") )
ggboxplot(df, x = "ID", y = "length",
          color = "dose", palette = "jco")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50) 

-------------------------------------------------------------------------------------------------------
#정규성 검사: Shapiro.test
df <- readRDS('/home/jurny0402/231116_nor,FTC_seurat/data/231212_cdr3length_distribution.rds')
df_adenoma <- subset(df, subset=ID=="a")
df_carcinoma <- subset(df, subset=ID=="c")
df_normal <- subset(df, subset=ID=="n")
df_ft5 <- subset(df_carcinoma, subset=values == "FT5_c")

df_a_length <- df_adenoma$length
df_c_length <- df_carcinoma$length
df_ft5_length <- df_ft5$length
df_n_length <- df_normal$length

write.table(df_a_length, "~/231116_nor,FTC_seurat/data/adenoma.txt", sep="\t", row.names = TRUE)
write.table(df_c_length, "~/231116_nor,FTC_seurat/data/carcinoma.txt", sep="\t", row.names = TRUE)
write.table(df_n_length, "~/231116_nor,FTC_seurat/data/normal.txt", sep="\t", row.names = TRUE)

result_a <- shapiro.test(df_a_length)
result_c <- shapiro.test(df_c_length)
result_n <- shapiro.test(df_n_length)
result_ft5 <- shapiro.test(df_ft5_length)

result_a;result_c;result_n;result_ft5
-------------------------------------------------------------------------------------------------------
#One-way ANOVA
library(ggplot2)
length(df_a_length)
length(df_c_length) 
length(df_n_length)
length(df_ft5_length)
  
merged_length <- data.frame(
  value= c(df_a_length, df_c_length, df_n_length),
  group = rep(c("adenoma","carcinoma","normal"), times = c(length(df_a_length), length(df_c_length), length(df_n_length)))
) 
  
result <- aov(value ~ group, data = merged_length)

summary(result)
-------------------------------------------------------------------------------------------------------  
#2 sample t test
ac.result <- t.test(df_a_length,df_c_length)
ac.result

an.result <- t.test(df_a_length,df_n_length)
an.result

cn.result <- t.test(df_c_length,df_n_length)
cn.result

data <- list(G1 = df_a_length, G2 = df_c_length, G3 = df_n_length)
result <- kruskal.test(data)
result
-------------------------------------------------------------------------------------------------------
#violin plot
head(df,4)
nrow(df)
for(i in 1:nrow(df)) {
  id <- df[i,3]
  if(id == "a"){
    df[i,3] <- c("adenoma")
  }
  else if(id == "c") {
    df[i,3] <-c("carcinoma")
  }
  else if(id == "n") {
    df[i,3] <- c("normal")
  }
}

library(ggplot2)
theme_set(
  theme_classic() +
    theme(legend.position = "top")
)  

e <- ggplot(df, aes(x = ID, y = length))

e + geom_violin(aes(fill = ID),trim = FALSE) + 
  stat_summary(
    fun.data = "mean_sdl",  fun.args = list(mult = 1), 
    geom = "pointrange", color = "black"
  )+
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))
  