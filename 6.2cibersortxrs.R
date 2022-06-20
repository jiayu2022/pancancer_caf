
getwd()

library(tidyverse)

library(patchwork)

library(dplyr)

library(limma)

library(ggplot2)

library(ggthemes)

library(ggpubr)

library(ggsci)

library(data.table)

rm(list = ls())

paad = read.csv(file = "3.1tcga1/1.paad/out/cellcafsub.csv")

paad$sampletype = rep("PAAD", nrow(paad))

table(paad$sampletype)

ovc = read.csv(file = "3.1tcga1/2.ovc/out/cellcafsub.csv")

ovc$sampletype = rep("OVC", nrow(ovc))

table(ovc$sampletype)

brca = read.csv(file = "3.1tcga1/3.brca/out/cellcafsub.csv")

brca$sampletype = rep("BRCA",nrow(brca))

table(brca$sampletype)

lihc = read.csv(file = "3.1tcga1/4.lihc/out/cellcafsub.csv")

lihc$sampletype = rep("LIHC", nrow(lihc))

table(lihc$sampletype)

coad = read.csv(file = "3.1tcga1/5.coad/out/cellcafsub.csv")

coad$sampletype = rep("COAD", nrow(coad))

table(coad$sampletype)

read = read.csv(file = "3.1tcga1/6.read/out/cellcafsub.csv")

read$sampletype = rep("READ", nrow(read))

table(read$sampletype)

prad = read.csv(file = "3.1tcga1/7.prad/out/cellcafsub.csv")

prad$sampletype = rep("PRAD", nrow(prad))

table(prad$sampletype)

luad = read.csv(file = "3.1tcga1/8.luad/out/cellcafsub.csv")

luad$sampletype = rep("LUAD", nrow(luad))

table(luad$sampletype)

lusc = read.csv(file = "3.1tcga1/9.lusc/out/cellcafsub.csv")

lusc$sampletype = rep("LUSC", nrow(lusc))

table(lusc$sampletype)

blca = read.csv(file = "3.2tcga2/1.blca/out/cellcafsub.csv")

blca$sampletype = rep("BLCA", nrow(blca))

table(blca$sampletype)

esca = read.csv(file = "3.2tcga2/2.esca/out/cellcafsub.csv")

esca$sampletype = rep("ESCA", nrow(esca))

table(esca$sampletype)

gbm = read.csv(file = "3.2tcga2/3.gbm/out/cellcafsub.csv")

gbm$sampletype = rep("GBM", nrow(gbm))

table(gbm$sampletype)

hnsc = read.csv(file = "3.2tcga2/4.hnsc/out/cellcafsub.csv")

hnsc$sampletype = rep("HNSC", nrow(hnsc))

table(hnsc$sampletype)

thca = read.csv(file = "3.2tcga2/5.thca/out/cellcafsub.csv")

thca$sampletype = rep("THCA", nrow(thca))

table(thca$sampletype)

stad = read.csv(file = "3.2tcga2/6.stad/out/cellcafsub.csv")

stad$sampletype = rep("STAD", nrow(stad))

table(stad$sampletype)

cesc = read.csv(file = "3.2tcga2/7.cesc/out/cellcafsub.csv")

cesc$sampletype = rep("CESC", nrow(cesc))

table(cesc$sampletype)

ucec = read.csv(file = "3.2tcga2/8.ucec/out/cellcafsub.csv")

ucec$sampletype = rep("UCEC", nrow(ucec))

table(ucec$sampletype)

colnames(paad)

colnames(ucec)

l = list(blca, brca, cesc, coad, esca,
         gbm, hnsc, lihc, luad, lusc,
         ovc, paad, prad, read, stad,
         thca, ucec)

panca = rbindlist(l, use.names=TRUE)

table(panca$sampletype)

#dir.create("figure/tcga")

#dir.create("figure/tcga/input")

dir.create("figure/tcga/out")

saveRDS(panca, file = "figure/tcga/input/pancact.rds")

colnames(panca)

table(panca$sample)

plot.info_all <- panca[, c(30, 22, 23)]

p1 = ggboxplot(plot.info_all,
               x = "sampletype",
               y = "FB",
               color = "sample",
               fill = "white", #只需要修改这里，讲fill = "CellType"改为fill = "group"
               outlier.size = 0.5,
               xlab = "",
               ylab = "Fraction of  Fibroblasts") +
  theme_classic() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5),
    axis.title.y = element_text(margin = margin(0, 0.35, 0, 0,'cm'))) +
  scale_color_d3() +
  stat_compare_means(aes(group = sample), label = "p.signif")

p1

ggsave("figure/tcga/out/fbsuball.png", plot = p1, width = 8, height = 4)

ggsave("figure/tcga/out/fbsuball.pdf", plot = p1, width = 8, height = 4)

colnames(panca)

plot.info_all2 <- panca[, c(30, 22, 24)]

p2 = ggboxplot(plot.info_all2,
               x = "sampletype",
               y = "COL11A1_FB.FB",
               color = "sample",
               fill = "white", #只需要修改这里，讲fill = "CellType"改为fill = "group"
               outlier.size = 0.5,
               xlab = "",
               ylab = "COL11A1+ FBs / Total FBs") +
  theme_classic() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5),
    axis.title.y = element_text(margin = margin(0, 0.35, 0, 0,'cm'))) +
  scale_color_d3() +
  stat_compare_means(aes(group = sample), label = "p.signif")

p2

ggsave("figure/tcga/out/col11a1fb_fb.png", plot = p2, width = 8, height = 4)

ggsave("figure/tcga/out/col11a1fb_fb.pdf", plot = p2, width = 8, height = 4)

colnames(panca)

plot.info_all3 <- panca[, c(30, 22, 16)]

p3 = ggboxplot(plot.info_all3,
               x = "sampletype",
               y = "COL11A1_FB",
               color = "sample",
               fill = "white", #只需要修改这里，讲fill = "CellType"改为fill = "group"
               outlier.size = 0.5,
               xlab = "",
               ylab = "Fraction of  COL11A1+ FBs") +
  theme_classic() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5),
    axis.title.y = element_text(margin = margin(0, 0.35, 0, 0,'cm'))) +
  scale_color_d3() +
  stat_compare_means(aes(group = sample), label = "p.signif")

p3

ggsave("figure/tcga/out/col11a1fb.png", plot = p3, width = 8, height = 4)

ggsave("figure/tcga/out/col11a1fb.pdf", plot = p3, width = 8, height = 4)



