
getwd()

library(Seurat)

library(tidyverse)

library(patchwork)

library(dplyr)

library("ggsci")

library(paletteer)

library(pheatmap) # 加载包

library(MySeuratWrappers)

library(Scillus)

rm(list = ls())

scRNA = readRDS(file = "keyobject/fbsubid.rds")

scRNA = readRDS(file = "input/keyobject/fbsubid.rds")

table(scRNA$fbsub)

table(scRNA$fbtype)

table(Idents(scRNA))

pal <- paletteer_d("ggsci::category20_d3")[c(5, 6, 2, 4, 3, 1)]

#dir.create("figure/marker2")

allgene = scRNA@assays$RNA@data@Dimnames[1]

str(allgene)

allgenes = allgene[[1]]

str(allgenes)

chemokines <- grep("^CXC|CCL|CCR|CX3|XCL|XCR", allgenes, value = T)

chemokines

p1 = VlnPlot(scRNA, features = chemokines, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
        theme(axis.text.x = element_text(angle = 90, size = rel(0.35))) +
        theme(axis.ticks.x = element_line(size = 0.35))

ggsave(p1, file = "figure/marker2/ckvioa.png", width = 20, height = 3.5)

ccl = grep("^CCL", allgenes, value = T)

ccl

str(ccl)

pal <- paletteer_d("ggsci::category20_d3")[c(5, 6, 2, 4, 3, 1)]

p1 = VlnPlot(scRNA, features = ccl, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
        theme(axis.text.x = element_text(angle = 90, size = rel(0.35))) +
        theme(axis.ticks.x = element_line(size = 0.35))

ggsave(p1, file = "figure/marker2/cclvioa.png", width = 10, height = 3.5)

cxcl = grep("^CXCL", allgenes, value = T)

cxcl

p1 = VlnPlot(scRNA, features = cxcl, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
        theme(axis.text.x = element_text(angle = 90, size = rel(0.35))) +
        theme(axis.ticks.x = element_line(size = 0.35))

ggsave(p1, file = "figure/marker2/cxclvioa.png", width = 6, height = 3.5)

il <- grep("^IL", allgenes, value = T)

il

p1 = VlnPlot(scRNA, features = il, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
        theme(axis.text.x = element_text(angle = 90, size = rel(0.35))) +
        theme(axis.ticks.x = element_line(size = 0.35))

ggsave(p1, file = "figure/marker2/ilvioa.png", width = 30, height = 3.5)

gf1 <- grep("^TGF|FGF|CTGF", allgenes, value = T)

gf2 <- grep("^IGF|HGF|EGF|VEGF|PDGF", allgenes, value = T)

gf1

gf2

p1 = VlnPlot(scRNA, features = gf1, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
        theme(axis.text.x = element_text(angle = 90, size = rel(0.35))) +
        theme(axis.ticks.x = element_line(size = 0.35))

ggsave(p1, file = "figure/marker2/gf1vioa.png", width = 30, height = 3.5)

p1 = VlnPlot(scRNA, features = gf2, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
        theme(axis.text.x = element_text(angle = 90, size = rel(0.35))) +
        theme(axis.ticks.x = element_line(size = 0.35))

ggsave(p1, file = "figure/marker2/gf2vioa.png", width = 30, height = 3.5)

col <- grep("^COL", allgenes, value = T)

col

p1 = VlnPlot(scRNA, features = col, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
        theme(axis.text.x = element_text(angle = 90, size = rel(0.35))) +
        theme(axis.ticks.x = element_line(size = 0.35))

ggsave(p1, file = "figure/marker2/colvioa.png", width = 30, height = 3.5)

mmp <- grep("^MMP", allgenes, value = T)

mmp

p1 = VlnPlot(scRNA, features = mmp, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
        theme(axis.text.x = element_text(angle = 90, size = rel(0.35))) +
        theme(axis.ticks.x = element_line(size = 0.35))

ggsave(p1, file = "figure/marker2/mmpvioa.png", width = 20, height = 3.5)

mk <- c("BAMBI", "CLDN1", "COL11A1", "CXCL14", "DPT", "RGS5")

p1 = VlnPlot(scRNA, features = mk, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
     theme(axis.text.x = element_blank(),
           axis.ticks.x = element_blank())  #不显示坐标刻度

p1 = VlnPlot(scRNA, features = mk, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "vertical", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())  #不显示坐标刻度

ggsave(p1, file = "figure/fig3e_mkvioa.png", width = 3.5, height = 4)

ggsave(p1, file = "figure/fig3e_mkvioa.pdf", width = 3.5, height = 4)

dir.create("figure/marker2022")

pal <- paletteer_d("ggsci::category20_d3")[c(5, 6, 2, 4, 3, 1)]

ECM = c("COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2",
        "COL5A1", "COL5A2", "COL6A1", "COL6A2", "COL6A3",
        "COL8A1", "COL8A2", "COL10A1", "COL11A1", "COL12A1",
        "COL14A1", "COL15A1", "COL16A1", "COL18A1",
        "MMP2", "MMP11", "MMP14", "MMP19", "MMP23B",
        "TIMP1", "TIMP2", "TIMP3",
        "FN1", "POSTN", "LOX")

p1 = VlnPlot(scRNA, features = ECM, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
     theme(axis.text.x = element_blank(),
           axis.ticks.x = element_blank())  #不显示坐标刻度

ggsave(p1, file = "figure/marker2022/ecmvio.png", width = 12, height = 3.5)

ggsave(p1, file = "figure/marker2022/ecmvio.pdf", width = 12, height = 3.5)

mk <- AverageExpression(scRNA,
                        group.by = "fbsub",
                        features = ECM,
                        return.seurat = FALSE)

mfdf = as.data.frame(mk[["RNA"]])

annotation_col = data.frame(Cluster = factor(c("BAMBI+", "CLDN1+", "COL11A1+", "CXCL14+", "DPT+", "RGS5+" )))

rownames(annotation_col)

colnames(mfdf)

rownames(annotation_col) <- colnames(mfdf)

head(annotation_col)

ann_colors = list(Cluster = c("DPT+" = "#2CA02CCC", "CLDN1+" = "#8C564BCC", "CXCL14+" = "#D62728CC",
                              "BAMBI+" = "#9467BDCC", "COL11A1+" = "#FF7F0ECC", "RGS5+" = "#1F77B4CC"))

head(ann_colors)

p <- pheatmap(mfdf,
              annotation_col = annotation_col,
              annotation_colors = ann_colors,
              annotation_legend = FALSE,
              annotation_names_col = FALSE,
              color = colorRampPalette(c("#0080FF", "white", "#FF8000"))(20),
              main = "ECM-associated",
              scale = "row",
              angle_col = 270, # 设置显示角度
              cellwidth = 20, cellheight = 20, # 设置热图方块宽度和高度
              border = "white",  # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F,
              show_rownames = T, #显示横、纵坐标id
              show_colnames = T,
              legend = T, # 显示图例
              fontsize_row = 10, # 分别设置横向和纵向字体大小
              fontsize_col = 10.5)

ggsave("figure/marker2022/ecm.png", p, width = 3, height = 10)

ggsave("figure/marker2022/ecm.pdf", p, width = 3, height = 10)

IM = c("IL18", "IL32", "IL33",
       "CCL2", "CCL8", "CCL11", "CCL13",
       "CXCL1", "CXCL2", "CXCL6", "CXCL12", "CXCL14",
       "CSF1")

p1 = VlnPlot(scRNA, features = IM, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
     theme(axis.text.x = element_blank(),
           axis.ticks.x = element_blank())  #不显示坐标刻度

ggsave(p1, file = "figure/marker2022/imvio.png", width = 8, height = 3.5)

ggsave(p1, file = "figure/marker2022/imvio.pdf", width = 8, height = 3.5)

mk <- AverageExpression(scRNA,
                        group.by = "fbsub",
                        features = IM,
                        return.seurat = FALSE)

mfdf = as.data.frame(mk[["RNA"]])

annotation_col = data.frame(Cluster = factor(c("BAMBI+", "CLDN1+", "COL11A1+", "CXCL14+", "DPT+", "RGS5+" )))

rownames(annotation_col)

colnames(mfdf)

rownames(annotation_col) <- colnames(mfdf)

head(annotation_col)

ann_colors = list(Cluster = c("DPT+" = "#2CA02CCC", "CLDN1+" = "#8C564BCC", "CXCL14+" = "#D62728CC",
                              "BAMBI+" = "#9467BDCC", "COL11A1+" = "#FF7F0ECC", "RGS5+" = "#1F77B4CC"))

head(ann_colors)

p <- pheatmap(mfdf,
              annotation_col = annotation_col,
              annotation_colors = ann_colors,
              annotation_legend = FALSE,
              annotation_names_col = FALSE,
              color = colorRampPalette(c("#0080FF", "white", "#FF8000"))(20),
              main = "Cytokine",
              scale = "row",
              angle_col = 270, # 设置显示角度
              cellwidth = 20, cellheight = 20, # 设置热图方块宽度和高度
              border = "white",  # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F,
              show_rownames = T, #显示横、纵坐标id
              show_colnames = T,
              legend = T, # 显示图例
              fontsize_row = 10, # 分别设置横向和纵向字体大小
              fontsize_col = 10.5)

ggsave("figure/marker2022/ck.png", p, width = 3, height = 5.5)

ggsave("figure/marker2022/ck.pdf", p, width = 3, height = 5.5)

GF = c("TGFB1", "TGFB3", "TGFBI", "TGFBR1", "TGFBR2", "TGFBR3",
       "FGF7", "FGFR1",
       "CTGF",
       "VEGFA", "VEGFB",
       "PDGFC", "PDGFD", "PDGFRA", "PDGFRB",
       "IGF1", "IGF2", "IGF1R", "IGF2R",
       "HGF",
       "EGFR")

p1 = VlnPlot(scRNA, features = GF, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
     theme(axis.text.x = element_blank(),
           axis.ticks.x = element_blank())  #不显示坐标刻度

ggsave(p1, file = "figure/marker2022/gfvio.png", width = 8, height = 3.5)

ggsave(p1, file = "figure/marker2022/gfvio.pdf", width = 8, height = 3.5)

mk <- AverageExpression(scRNA,
                        group.by = "fbsub",
                        features = GF,
                        return.seurat = FALSE)

mfdf = as.data.frame(mk[["RNA"]])

p <- pheatmap(mfdf,
              annotation_col = annotation_col,
              annotation_colors = ann_colors,
              annotation_legend = FALSE,
              annotation_names_col = FALSE,
              color = colorRampPalette(c("#0080FF", "white", "#FF8000"))(20),
              main = "Growth Factor",
              scale = "row",
              angle_col = 270, # 设置显示角度
              cellwidth = 20, cellheight = 20, # 设置热图方块宽度和高度
              border = "white",  # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F,
              show_rownames = T, #显示横、纵坐标id
              show_colnames = T,
              legend = T, # 显示图例
              fontsize_row = 10, # 分别设置横向和纵向字体大小
              fontsize_col = 10.5)

ggsave("figure/marker2022/gf.png", p, width = 3, height = 7.5)

ggsave("figure/marker2022/gf.pdf", p, width = 3, height = 7.5)

markers <- c("CD248",
             "PDGFRA",
             "FAP", "PDPN",
             "LUM", "COL6A3", "THY1", "DCN",
             "COL3A1", "COL1A1", "COL1A2")

p1 = VlnPlot(scRNA, features = markers, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
     theme(axis.text.x = element_blank(),
           axis.ticks.x = element_blank())  #不显示坐标刻度

ggsave(p1, file = "figure/marker2022/fbmkvio.png", width = 5, height = 3.5)

ggsave(p1, file = "figure/marker2022/fbmkvio.pdf", width = 5, height = 3.5)

markers <- c("IGF1", "IGF2", "HGF", "WNT5A", "HBEGF",
             "VEGFA", "VEGFB", "PGF",
             "FGF2", "FGF7",
             "TGFB1", "TGFB3",
             "PDGFA", "PDGFC", "PDGFD",
             "GAS6", "SPP1")

p1 = VlnPlot(scRNA, features = markers, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
     theme(axis.text.x = element_blank(),
           axis.ticks.x = element_blank())  #不显示坐标刻度

ggsave(p1, file = "figure/marker2022/fblrvio.png", width = 8, height = 3.5)

ggsave(p1, file = "figure/marker2022/fblrvio.pdf", width = 8, height = 3.5)

markers <- c("IGF1", "IGF2", "HGF", "WNT5A",
             "VEGFA", "VEGFB", "PGF",
             "FGF7",
             "TGFB1", "TGFB3",
             "PDGFC", "PDGFD",
             "GAS6")

p1 = VlnPlot(scRNA, features = markers, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
     theme(axis.text.x = element_blank(),
           axis.ticks.x = element_blank())  #不显示坐标刻度

ggsave(p1, file = "figure/marker2022/fblrvioa.png", width = 6.5, height = 3.5)

ggsave(p1, file = "figure/marker2022/fblrvioa.pdf", width = 6.5, height = 3.5)

markers <- c("FAP", "LRRC15", "ITGA11",
             "SDC1", "SPHK1", "ADAM12",
             "TNFRSF12A", "IFITM1", "ITGB5", "ITGB1", "THY1")

p1 = VlnPlot(scRNA, features = markers, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
     theme(axis.text.x = element_blank(),
           axis.ticks.x = element_blank())  #不显示坐标刻度

ggsave(p1, file = "figure/marker2022/fbmembvioa.png", width = 5, height = 3.5)

ggsave(p1, file = "figure/marker2022/fbmembvioa.pdf", width = 5, height = 3.5)

markers <- c("POSTN", "COL11A1", "COL5A2")

p1 = VlnPlot(scRNA, features = markers, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "vertical",
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
     theme(axis.text.y = element_blank(),
           axis.ticks.y = element_blank())  #不显示坐标刻度

ggsave(p1, file = "figure/marker2022/fbsubkmvioa.png", width = 5, height = 3.5)

ggsave(p1, file = "figure/marker2022/fbsubkmvioa.pdf", width = 5, height = 3.5)

markers <- c("POSTN", "COL11A1", "COL5A2")

pal <- paletteer_d("ggsci::category20_d3")[c(2, 1)]

p1 = VlnPlot(scRNA, features = markers, cols = pal, group.by = "fbtype",
             stacked = T, pt.size = 0,
             direction = "vertical",
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())  #不显示坐标刻度

p1

ggsave(p1, file = "figure/fig6b.png", width = 2.8, height = 3.5)

ggsave(p1, file = "figure/fig6b.pdf", width = 2.8, height = 3.5)



