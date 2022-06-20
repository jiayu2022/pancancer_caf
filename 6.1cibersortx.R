
setwd("/home/zhangjiayu/project/hucafnf")

getwd()

library(Seurat)

library(tidyverse)

library(patchwork)

library(dplyr)

rm(list = ls())

scRNA = readRDS(file = "keyobject/fbsubid.rds")

table(Idents(scRNA))

str(scRNA@meta.data)

table(scRNA$fbtype)

Cells.sub <- subset(scRNA@meta.data, fbtype == "CAF")

caf <- subset(scRNA, cells = row.names(Cells.sub))

table(caf$fbtype)

saveRDS(caf, file = "keyobject/cafsubid.rds")

