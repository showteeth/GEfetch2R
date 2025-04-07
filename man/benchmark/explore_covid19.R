# library
library(Seurat)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(GEfetch2R)

# prepare the data
# ucsc.cb.samples = ShowCBDatasets(lazy = TRUE, json.folder = "/Volumes/soyabean/GEfetch2R/cell_browser/json", update = TRUE)
# saveRDS(ucsc.cb.samples, file = "ucsc.cb.samples.rds")
# ucsc.cb.samples.rds is available via https://github.com/showteeth/GEfetch2R/tree/main/man/benchmark
ucsc.cb.samples = readRDS('ucsc.cb.samples.rds')
# extract information of selected dataset
covidT.sample.df <- ExtractCBDatasets(all.samples.df = ucsc.cb.samples, collection = "COVID-19 Immunological Response",
                                      sub.collection = "T-cells")
# extract cell type compositional
covidT.sample.ct <- ExtractCBComposition(json.folder = "/Volumes/soyabean/GEfetch2R/cell_browser/json", sample.df = covidT.sample.df)
# download dataset and load to Seurat
covidT.sample.seu <- ParseCBDatasets(sample.df = covidT.sample.df)
# get valid data
covidT.sample.seu.meta = covidT.sample.seu@meta.data
covidT.sample.seu.meta = covidT.sample.seu.meta[!is.na(covidT.sample.seu.meta$UMAP_1), ]

# Fig. 2A
covidT.cluster.plot = ggplot() +
  geom_point(data = covidT.sample.seu.meta, aes(x = UMAP_1, y=UMAP_2, color = FinalCellTypeTcell)) +
  scale_color_manual(values = clu.colors) +
  theme_classic() +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = .6))
ggsave(filename = "covidT_cluster_plot.pdf", plot = covidT.cluster.plot, height = 5, width = 7)

# proportion of each T cell subtype 
covidT.sample.seu.meta.stat = covidT.sample.seu.meta %>% 
  dplyr::group_by(orig.ident, FinalCellTypeTcell) %>% 
  dplyr::summarise(Num = n()) %>% 
  dplyr::group_by(orig.ident) %>% 
  dplyr::mutate(Freq = 100 * Num/sum(Num))
# Fig. 2B
t.prop.plot = ggplot(covidT.sample.seu.meta.stat, aes(x=orig.ident, y=Freq, fill = FinalCellTypeTcell)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values = clu.colors) +
  theme_classic() +
  theme(axis.title.x = element_blank()) + 
  scale_y_continuous(limits = c(0, 100), expand = expansion(add = 1)) + 
  scale_x_discrete(expand = expansion(add = 0.5))
ggsave(filename = "covidT_t_proportion.pdf", plot = t.prop.plot, height = 5, width = 6)

# summarize the data
covidT.sample.ct.stat = covidT.sample.seu.meta %>%
  dplyr::group_by(PatientID, orig.ident, FinalCellTypeTcell) %>%
  dplyr::summarise(Num = n()) %>%
  dplyr::mutate(Freq = Num / sum(Num)) %>%
  as.data.frame()
# perform Mann–Whitney U-test
covidT.sample.ct.stat$FinalCellTypeTcell = 
  factor(covidT.sample.ct.stat$FinalCellTypeTcell, levels = c("CD4+ Naive", "CD4+ Memory", "CD4+ Effector Memory", "CD4+ Effector-GZMK", "CD4+ Effector-GNLY", "Treg",
                                                              "CD8+ Naive", "CD8+ Effector-GZMK", "CD8+ Effector-GNLY", "NKT Naive", "NKT CD56", "NKT CD160"))
covidT.sample.ct.stat.test = covidT.sample.ct.stat %>%
  dplyr::group_by(FinalCellTypeTcell) %>%
  rstatix::wilcox_test(Freq ~ orig.ident) %>% 
  dplyr::filter(p < 0.05) %>%
  rstatix::add_xy_position(x = "FinalCellTypeTcell", step.increase = 0.02)
# Fig. 2C
clu.group.compositional = ggplot(covidT.sample.ct.stat, aes(x=FinalCellTypeTcell, y=Freq, fill = orig.ident)) +
  geom_bar(stat="summary", fun = mean, position='dodge') +
  geom_point(aes(x = FinalCellTypeTcell), position=
               position_jitterdodge(dodge.width=0.9), color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.4, 
               position = position_dodge(0.9)) + 
  scale_fill_manual(values = sample.colors) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggpubr::stat_pvalue_manual(covidT.sample.ct.stat.test, label = "p", tip.length = 0, inherit.aes = F) + 
  scale_y_continuous(limits = c(0, 0.55), expand = expansion(add = 0.01))
ggsave(filename = "covidT_cluster_group_compositional.pdf", plot = clu.group.compositional, height = 5, width = 10)

# perform unpaired Dunn’s test
covidT.sample.seu.meta$FinalCellTypeTcell = 
  factor(covidT.sample.seu.meta$FinalCellTypeTcell, levels = c("CD4+ Naive", "CD4+ Memory", "CD4+ Effector Memory", "CD4+ Effector-GZMK", "CD4+ Effector-GNLY", "Treg",
                                                               "CD8+ Naive", "CD8+ Effector-GZMK", "CD8+ Effector-GNLY", "NKT Naive", "NKT CD56", "NKT CD160"))
covidT.sample.seu.meta.apop.test = covidT.sample.seu.meta %>%
  group_by(FinalCellTypeTcell) %>%
  rstatix::dunn_test(Apoptosis_score ~ orig.ident, p.adjust.method = "bonferroni") %>% 
  dplyr::filter(p.adj < 0.01) %>%
  rstatix::add_xy_position(x = "FinalCellTypeTcell", step.increase = 0.08)
# Fig. 2D
apop.score.plot = ggplot() +
  geom_boxplot(data = covidT.sample.seu.meta,
               aes(x = FinalCellTypeTcell, y=Apoptosis_score,
                   fill = orig.ident), outlier.shape = NA) +
  scale_fill_manual(values = sample.colors) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggpubr::stat_pvalue_manual(covidT.sample.seu.meta.apop.test, label = "p.adj.signif", tip.length = 0) + 
  scale_y_continuous(limits = c(0, 0.25), expand = expansion(add = 0.01))
ggsave(filename = "covidT_apop_score.pdf", plot = apop.score.plot, height = 5, width = 10)




