
getwd()

library(tidyverse)

library(patchwork)

library(dplyr)

library(ggplot2)

library(ggthemes)

library(ggsci)

library(data.table)

library(scales)

rm(list = ls())

paad = read.csv(file = "3.1tcga1/1.paad/coxkm/sctcoxa.csv")

paad$sampletype = rep("PAAD", nrow(paad))

table(paad$sampletype)

ovc = read.csv(file = "3.1tcga1/2.ovc/coxkm/sctcoxa.csv")

ovc$sampletype = rep("OVC", nrow(ovc))

table(ovc$sampletype)

brca = read.csv(file = "3.1tcga1/3.brca/coxkm/sctcoxa.csv")

brca$sampletype = rep("BRCA",nrow(brca))

table(brca$sampletype)

lihc = read.csv(file = "3.1tcga1/4.lihc/coxkm/sctcoxa.csv")

lihc$sampletype = rep("LIHC", nrow(lihc))

table(lihc$sampletype)

coad = read.csv(file = "3.1tcga1/5.coad/coxkm/sctcoxa.csv")

coad$sampletype = rep("COAD", nrow(coad))

table(coad$sampletype)

read = read.csv(file = "3.1tcga1/6.read/coxkm/sctcoxa.csv")

read$sampletype = rep("READ", nrow(read))

table(read$sampletype)

prad = read.csv(file = "3.1tcga1/7.prad/coxkm/sctcoxa.csv")

prad$sampletype = rep("PRAD", nrow(prad))

table(prad$sampletype)

luad = read.csv(file = "3.1tcga1/8.luad/coxkm/sctcoxa.csv")

luad$sampletype = rep("LUAD", nrow(luad))

table(luad$sampletype)

lusc = read.csv(file = "3.1tcga1/9.lusc/coxkm/sctcoxa.csv")

lusc$sampletype = rep("LUSC", nrow(lusc))

table(lusc$sampletype)

blca = read.csv(file = "3.2tcga2/1.blca/coxkm/sctcoxa.csv")

blca$sampletype = rep("BLCA", nrow(blca))

table(blca$sampletype)

esca = read.csv(file = "3.2tcga2/2.esca/coxkm/sctcoxa.csv")

esca$sampletype = rep("ESCA", nrow(esca))

table(esca$sampletype)

gbm = read.csv(file = "3.2tcga2/3.gbm/coxkm/sctcoxa.csv")

gbm$sampletype = rep("GBM", nrow(gbm))

table(gbm$sampletype)

hnsc = read.csv(file = "3.2tcga2/4.hnsc/coxkm/sctcoxa.csv")

hnsc$sampletype = rep("HNSC", nrow(hnsc))

table(hnsc$sampletype)

thca = read.csv(file = "3.2tcga2/5.thca/coxkm/sctcoxa.csv")

thca$sampletype = rep("THCA", nrow(thca))

table(thca$sampletype)

stad = read.csv(file = "3.2tcga2/6.stad/coxkm/sctcoxa.csv")

stad$sampletype = rep("STAD", nrow(stad))

table(stad$sampletype)

cesc = read.csv(file = "3.2tcga2/7.cesc/coxkm/sctcoxa.csv")

cesc$sampletype = rep("CESC", nrow(cesc))

table(cesc$sampletype)

ucec = read.csv(file = "3.2tcga2/8.ucec/coxkm/sctcoxa.csv")

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

dir.create("figure/tcga/coxkm")

saveRDS(panca, file = "figure/tcga/coxkm/pancact.rds")

colnames(panca)

table(panca$sample)

panca$pstar <- ifelse(panca$pvalue < 0.05,
                      ifelse(panca$pvalue < 0.01,
                             ifelse(panca$pvalue < 0.001,"***", "**"), "*"), "")

write.csv(panca, file = "figure/tcga/coxkm/pancact.csv")

panca = read.csv(file = "figure/tcga/coxkm/pancacta.csv")

saveRDS(panca, file = "figure/tcga/coxkm/pancact.rds")

panca = readRDS(file = "figure/tcga/coxkm/pancact.rds")

p1 = ggplot(panca, aes(sampletype, id)) +
        geom_tile(aes(fill = HR), colour = "white", size = 0.1) +
        scale_fill_gradient2(low = "#0080FF", mid = "white", high = "red", midpoint = 1) +
        geom_text(aes(label = pstar), col ="black", size = 5) +
        theme_minimal() +
        theme(axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              axis.text.y = element_text(size = 8))

p1

ggsave(p1, file = "figure/tcga/coxkm/cox.png", width = 6.5, height = 2.8)

ggsave(p1, file = "figure/tcga/coxkm/cox.pdf", width = 6.5, height = 2.8)



