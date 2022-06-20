
getwd()

library(ggplot2)

library(ggthemes)

library(dplyr)

library(data.table)

library(corrplot)

rm(list = ls())

paad = readRDS(file = "3.1tcga1/1.paad/immu/correlation.rds")

paad$cancertype = rep("PAAD", nrow(paad))

table(paad$cancertype)

ovc = readRDS(file = "3.1tcga1/2.ovc/immu/correlation.rds")

ovc$cancertype = rep("OVC", nrow(ovc))

table(ovc$cancertype)

brca = readRDS(file = "3.1tcga1/3.brca/immu/correlation.rds")

brca$cancertype = rep("BRCA", nrow(brca))

table(brca$cancertype)

lihc = readRDS(file = "3.1tcga1/4.lihc/immu/correlation.rds")

lihc$cancertype = rep("LIHC", nrow(lihc))

table(lihc$cancertype)

coad = readRDS(file = "3.1tcga1/5.coad/immu/correlation.rds")

coad$cancertype = rep("COAD", nrow(coad))

table(coad$cancertype)

read = readRDS(file = "3.1tcga1/6.read/immu/correlation.rds")

read$cancertype = rep("READ", nrow(read))

table(read$cancertype)

prad = readRDS(file = "3.1tcga1/7.prad/immu/correlation.rds")

prad$cancertype = rep("PRAD", nrow(prad))

table(prad$cancertype)

luad = readRDS(file = "3.1tcga1/8.luad/immu/correlation.rds")

luad$cancertype = rep("LUAD", nrow(luad))

table(luad$cancertype)

lusc = readRDS(file = "3.1tcga1/9.lusc/immu/correlation.rds")

lusc$cancertype = rep("LUSC", nrow(lusc))

table(lusc$cancertype)

blca = readRDS(file = "3.2tcga2/1.blca/immu/correlation.rds")

blca$cancertype = rep("BLCA", nrow(blca))

table(blca$cancertype)

esca = readRDS(file = "3.2tcga2/2.esca/immu/correlation.rds")

esca$cancertype = rep("ESCA", nrow(esca))

table(esca$cancertype)

gbm = readRDS(file = "3.2tcga2/3.gbm/immu/correlation.rds")

gbm$cancertype = rep("GBM", nrow(gbm))

table(gbm$cancertype)

hnsc = readRDS(file = "3.2tcga2/4.hnsc/immu/correlation.rds")

hnsc$cancertype = rep("HNSC", nrow(hnsc))

table(hnsc$cancertype)

thca = readRDS(file = "3.2tcga2/5.thca/immu/correlation.rds")

thca$cancertype = rep("THCA", nrow(thca))

table(thca$cancertype)

stad = readRDS(file = "3.2tcga2/6.stad/immu/correlation.rds")

stad$cancertype = rep("STAD", nrow(stad))

table(stad$cancertype)

cesc = readRDS(file = "3.2tcga2/7.cesc/immu/correlation.rds")

cesc$cancertype = rep("CESC", nrow(cesc))

table(cesc$cancertype)

ucec = readRDS(file = "3.2tcga2/8.ucec/immu/correlation.rds")

ucec$cancertype = rep("UCEC", nrow(ucec))

table(ucec$cancertype)

colnames(paad)

colnames(ucec)

l = list(blca, brca, cesc, coad, esca,
         gbm, hnsc, lihc, luad, lusc,
         ovc, paad, prad, read, stad,
         thca, ucec)

panca = rbindlist(l, use.names=TRUE)

table(panca$cancertype)

dir.create("figure/tcga/cellfrac")

colnames(panca)

panca$pstar <- ifelse(panca$p.value < 0.05,
                     ifelse(panca$p.value < 0.01,
                            ifelse(panca$p.value < 0.001,"***", "**"), "*"), "")

write.csv(panca, file = "figure/tcga/cellfrac/pancact.csv")

panca = read.csv(file = "figure/tcga/cellfrac/pancacta.csv")

table(panca$ieep_cells)

panca = filter(panca, ieep_cells != "Epithelial_cells")

table(panca$ieep_cells)

saveRDS(panca, file = "figure/tcga/cellfrac/pancact.rds")

pancatb = tibble::as_tibble(panca)

table(panca$factori)

col11a1 = filter(pancatb, factori == "COL11A1_FB")

p1 = ggplot(col11a1, aes(cancertype, ieep_cells)) +
  geom_tile(aes(fill = cor), colour = "white", size = 0.1) +
  scale_fill_gradient2(low = "#0080FF", mid = "white", high = "red") +
  geom_text(aes(label = pstar), col ="black", size = 5) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8))+
  labs(title = "COL11A1+ FB", fill =paste0("Correlation")) +
  theme(plot.title = element_text(hjust = 0.5, size = 10.5))

p1

ggsave(p1, file = "figure/tcga/cellfrac/col11a1fbcor.png", width = 6.5, height = 3.2)

ggsave(p1, file = "figure/tcga/cellfrac/col11a1fbcor.pdf", width = 6.5, height = 3.2)

cxcl14 = filter(pancatb, factori == "CXCL14_FB")

p1 = ggplot(cxcl14, aes(cancertype, ieep_cells)) +
  geom_tile(aes(fill = cor), colour = "white", size = 0.1) +
  scale_fill_gradient2(low = "#0080FF", mid = "white", high = "red") +
  geom_text(aes(label = pstar), col ="black", size = 5) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8))+
  labs(title = "CXCL14+ FB", fill =paste0("Correlation")) +
  theme(plot.title = element_text(hjust = 0.5, size = 10.5))

p1

ggsave(p1, file = "figure/tcga/cellfrac/cxcl14fbcor.png", width = 6.5, height = 3.2)

ggsave(p1, file = "figure/tcga/cellfrac/cxcl14fbcor.pdf", width = 6.5, height = 3.2)

rgs5 = filter(pancatb, factori == "RGS5_FB")

p1 = ggplot(rgs5, aes(cancertype, ieep_cells)) +
  geom_tile(aes(fill = cor), colour = "white", size = 0.1) +
  scale_fill_gradient2(low = "#0080FF", mid = "white", high = "red") +
  geom_text(aes(label = pstar), col ="black", size = 5) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8))+
  labs(title = "RGS5+ FB", fill =paste0("Correlation")) +
  theme(plot.title = element_text(hjust = 0.5, size = 10.5))

p1

ggsave(p1, file = "figure/tcga/cellfrac/RGS5fbcor.png", width = 6.5, height = 3.2)

ggsave(p1, file = "figure/tcga/cellfrac/RGS5fbcor.pdf", width = 6.5, height = 3.2)

table(pancatb$factori)

dpt = filter(pancatb, factori == "DPT_FB")

p1 = ggplot(dpt, aes(cancertype, ieep_cells)) +
  geom_tile(aes(fill = cor), colour = "white", size = 0.1) +
  scale_fill_gradient2(low = "#0080FF", mid = "white", high = "red") +
  geom_text(aes(label = pstar), col ="black", size = 5) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8))+
  labs(title = "DPT+ FB", fill =paste0("Correlation")) +
  theme(plot.title = element_text(hjust = 0.5, size = 10.5))

p1

ggsave(p1, file = "figure/tcga/cellfrac/DPTfbcor.png", width = 6.5, height = 3.2)

ggsave(p1, file = "figure/tcga/cellfrac/DPTfbcor.pdf", width = 6.5, height = 3.2)

cldn1 = filter(pancatb, factori == "CLDN1_FB")

p1 = ggplot(cldn1, aes(cancertype, ieep_cells)) +
  geom_tile(aes(fill = cor), colour = "white", size = 0.1) +
  scale_fill_gradient2(low = "#0080FF", mid = "white", high = "red") +
  geom_text(aes(label = pstar), col ="black", size = 5) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8))+
  labs(title = "CLDN1+ FB", fill =paste0("Correlation")) +
  theme(plot.title = element_text(hjust = 0.5, size = 10.5))

p1

ggsave(p1, file = "figure/tcga/cellfrac/cldn1fbcor.png", width = 6.5, height = 3.2)

ggsave(p1, file = "figure/tcga/cellfrac/cldn1fbcor.pdf", width = 6.5, height = 3.2)

bambi = filter(pancatb, factori == "BAMBI_FB")

p1 = ggplot(bambi, aes(cancertype, ieep_cells)) +
  geom_tile(aes(fill = cor), colour = "white", size = 0.1) +
  scale_fill_gradient2(low = "#0080FF", mid = "white", high = "red") +
  geom_text(aes(label = pstar), col ="black", size = 5) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8))+
  labs(title = "BAMBI+ FB", fill =paste0("Correlation")) +
  theme(plot.title = element_text(hjust = 0.5, size = 10.5))

p1

ggsave(p1, file = "figure/tcga/cellfrac/BAMBIfbcor.png", width = 6.5, height = 3.2)

ggsave(p1, file = "figure/tcga/cellfrac/BAMBIfbcor.pdf", width = 6.5, height = 3.2)
