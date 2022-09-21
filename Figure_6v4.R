# whole_exome sequencing (actually for figure 6)
# only use exome sequencing index 
# Figure 6: 

# TCGA sig
#nature re-fitting (tissue color-coded)
# TCGA Val vs AA enrichment (tissue color-coded)
# TCGA Hypermutators
# TCGA POLE status
# HMF WES refitting
# HMF WES AA enrichment
# Confusion table HMF data (replaced in this version by a venn-diagram)

setwd("E:/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/")
source("E:/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Manuscript/Figures/Figure_4/Figure_4_v3.R")


# get samples with somatic POLE mutations 
POLE_muts = list.files("New_HMF/POLE_hotspot_muts/", full.names = T)
names(POLE_muts) = gsub("_.*", "", basename(POLE_muts))
POLE_changes = sapply(POLE_muts, function(x) strsplit(read.table(x)[,8], "\\|")[[1]][11])
metadata$POLE_somatic = "WT"
metadata$POLE_somatic[match(names(POLE_changes), metadata$sampleId)] = POLE_changes
metadata = metadata %>% mutate(POLE_somatic = replace(POLE_somatic,!grepl("Val411Leu|Pro286Arg|Ser459Phe", POLE_somatic), "WT"))
metadata$POLE_somatic = ifelse(metadata$POLE_somatic != "WT", "MUT", "WT")
metadata$log10muts = log10(metadata$n_sbs)

meta_exome = metadata[exome_select,]
meta_exome$SBS88_WES = 0
meta_exome$ID18_WES = 0
meta_exome[colnames(sbs_contri_exome), "SBS88_WES"] = sbs_contri_exome["SBS88",]
meta_exome[colnames(id_contri_exome), "ID18_WES"] = id_contri_exome[ "ID18",]
meta_exome$WES_sig = ifelse(meta_exome$SBS88_WES > 0.05 & meta_exome$ID18_WES > 0.05, "sig+", "sig-")
meta_exome$MS_WES = ifelse(meta_exome$log_p_all_exome > 3 & meta_exome$AA_exome > 0.26, "motif+","motif-")
meta_exome$WES_motif_sig = paste0(meta_exome$MS_WES, "/", meta_exome$WES_sig)

motifpos = meta_exome$sampleId[meta_exome$Motif_selection == "motif+"]
sigpos = meta_exome$sampleId[meta_exome$Sig_classification == "sig+"]
sigpos_WES = meta_exome$sampleId[meta_exome$WES_sig == "sig+"]
motifpos_WES = meta_exome$sampleId[meta_exome$MS_WES == "motif+"]

meta_exome[meta_exome$Sig_classification == "sig+", "primaryTumorLocation"] %>% table()
meta_exome[meta_exome$Motif_selection == "motif+", "primaryTumorLocation"] %>% table()



sigs_exome_HMF = ggplot(meta_exome, aes(x = ID18_WES, y = SBS88_WES, color = Motif_selection, shape = POLE_somatic)) + 
  geom_point(alpha = 0.7) + 
 # scale_color_manual(values = c("lightblue", "#d95f02", 'black', 'green'), name = "WGS motif selection") + 
  scale_shape_discrete(name = "Pks classification") + 
  scale_shape_manual(values = c(17, 16))+ 
  xlab("ID18 relative contribution") + ylab("SBS 88 relative contribution") + 
  theme_classic() + ggtitle("Signature refitting\n HMF exonic regions") + 
  theme(plot.title = element_text(hjust = 0.5))

sigs_exome_HMF
ggsave("../Manuscript/Figures/Figure_6/sig_exome.pdf", sigs_exome_HMF, width = 7, height = 4)

p_value_AA_exome_HMF = ggplot(meta_exome, aes(x = AA_exome, y = log_p_all_exome, color = Motif_selection, shape = POLE_somatic)) +
  geom_point(alpha = 0.7) + 
  geom_hline(yintercept = 3) + geom_vline(xintercept = 0.26) + 
#  scale_color_manual(values = c("lightblue", "#d95f02", 'black', 'green'), name = "WGS motif selection") + 
  scale_shape_manual(values = c(17, 16))+ 
  xlab("fraction mutations with AA at\n-3-4 position at pks-sites") + ylab("-log10 p-value") + 
  theme_classic() + ggtitle("motif selection HMF exonic regions") + 
  theme(plot.title = element_text(hjust = 0.5))

p_value_AA_exome_HMF
ggsave("../Manuscript/Figures/Figure_6/p_value_AA_exome.pdf", p_value_AA_exome_HMF, width = 7, height = 4)


unique(meta_exome$rec_tissues )
p_value_AA_exome_HMF_tissue = ggplot(meta_exome, aes(x = AA_exome, y = log_p_all_exome, color = rec_tissues, shape = POLE_somatic)) +
  geom_point(alpha = 0.7) + 
  geom_hline(yintercept = 3) + geom_vline(xintercept = 0.26) + 
  scale_color_manual(values = c("lightblue", "#d95f02", "#1b9e77",  "#bdbdbd",
                                "#7570b3", "black", "maroon", 'pink', 'violet'), name = "Primary cancer origin") + 
  scale_shape_manual(values = c(17, 16))+ 
  xlab("fraction mutations with AA at\n-3-4 position at pks-sites") + ylab("-log10 p-value") + 
  theme_classic() + ggtitle("motif selection HMF exonic regions") + 
  theme(plot.title = element_text(hjust = 0.5))
p_value_AA_exome_HMF_tissue
ggsave("../Manuscript/Figures/Figure_6/p_value_AA_exome_tissue.pdf", p_value_AA_exome_HMF_tissue, width = 7, height = 4)


# # include the plots by tissue 
# HMF_h_pval = ggplot(meta_exome, aes(x = rec_tissues, y = log_p_all_exome, color = Motif_Sig, shape = POLE_somatic)) +
#   scale_color_manual(values = c("lightblue", "#d95f02", 'black', 'green'), name = "WGS motif selection") + 
#   geom_jitter(alpha = 0.7, size = 1.7, width = 0.3) + 
#   geom_hline(yintercept = 3) + geom_vline(xintercept = 0.26) + 
#   scale_fill_manual(values = c("lightblue", "#d95f02", 'black', 'green'), name = "WGS motif selection") + 
#   scale_shape_manual(values = c(17, 16))+ 
#   theme_classic() +  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + xlab("") + 
#   ylab("-log10 p-value")
# 
# HMF_h_AA = ggplot(meta_exome, aes(x = rec_tissues, y = AA, color = Motif_Sig, shape = POLE_somatic )) +
#   geom_jitter(alpha = 0.7, size = 1.7, width = 0.3) + 
#   geom_vline(xintercept = 0.26) + 
#   scale_fill_manual(values = c("lightblue", "#d95f02", 'black', 'green'), name = "WGS motif selection") + 
#   scale_shape_manual(values = c(17, 16))+ 
#   theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + xlab("") + 
#   ylab("AA fraction")


# venn diagram comparing the signature identification methods: 
venn_list = list(`WGS motif` = motifpos, `WGS Signature` = sigpos, `WES motif` = motifpos_WES, `WES signature` = sigpos_WES)
venn = VennDiagram::venn.diagram(x = venn_list,filename = NULL,
                                 col= c('#8c510a','#d8b365','#5ab4ac','#01665e'), 
                                 fill = c(alpha("#8c510a", 0.3), alpha("#d8b365", 0.3), 
                                          alpha("#5ab4ac", 0.3), alpha("#01665e", 0.3))) %>% 
  grid::grobTree() %>% as_ggplot()

ggsave("../Manuscript/Figures/Figure_6/Vennplot.pdf", venn, width = 6, height = 6)


# Figure 6 updated with TCGA data Henry Wood 
meta_TCGA = read.delim(file = "../Other_datasets/TCGA/TCGA_metadata_noID.tsv")
meta_TCGA$Motif_selection = ifelse(meta_TCGA$pval_dinucs > 3 & meta_TCGA$AA > 0.26, "motif+", "motif-")
meta_TCGA$POLE_somatic = factor(meta_TCGA$POLE_somatic, levels = c("WT", "MUT"))

sigs_exome_ms = ggplot(meta_TCGA, aes(x = ID18, y = SBS88, color = cohort, shape = POLE_somatic)) +
  geom_point(alpha = 0.7) + 
  scale_shape_manual(values = c(16, 17))+ 
  geom_hline(yintercept = 0.05) + geom_vline(xintercept = 0.05) + 
  scale_color_manual(values =  c("#1b9e77",  "#bdbdbd","#7570b3", "black")) + 
  xlab("ID18 relative contribution") + ylab("SBS 88 relative contribution") + 
  theme_classic() + ggtitle("TCGA signature refitting\n exonic regions") + 
  theme(plot.title = element_text(hjust = 0.5))
sigs_exome_ms
ggsave("../Manuscript/Figures/Figure_6/TCGA_sig_exome.pdf", sigs_exome_ms, width = 7, height = 4)

p_value_AA_exome_ms = ggplot(meta_TCGA, aes(x = AA, y = pval_dinucs, color =  cohort, shape = POLE_somatic)) + 
  geom_point(alpha = 0.7) + 
  geom_hline(yintercept = 3) + geom_vline(xintercept = 0.26) + 
  scale_shape_manual(values = c(16, 17)) +
  scale_color_manual(values =  c("#1b9e77",  "#bdbdbd","#7570b3", "black")) + 
  xlab("fraction mutations with AA at\n-3-4 position at pks-sites") +
  ylab("-log10 p-value") + 
  theme_classic() + ggtitle("TCGA_motif selection exonic regions") + 
  theme(plot.title = element_text(hjust = 0.5))
p_value_AA_exome_ms
ggsave("../Manuscript/Figures/Figure_6/TCGA_p_value_AA_exome.pdf", p_value_AA_exome_ms, width = 7, height = 4)

TCGA_h_pval = ggplot(meta_TCGA, aes(x = cohort, y =pval_dinucs , color =  Motif_selection, shape = POLE_somatic)) +
  geom_jitter(alpha = 0.7, size = 1.7) + 
  scale_shape_manual(values = c(16, 17)) + 
  xlab("TCGA cohort") + ylab("-log10 p-value") + 
  theme_classic() + ggtitle("TCGA_motif selection exonic regions") + 
  theme(plot.title = element_text(hjust = 0.5))

TCGA_h_aa = ggplot(meta_TCGA, aes(x = cohort, y =AA , color =  Motif_selection, shape = POLE_somatic)) +
  geom_jitter(alpha = 0.7, size = 1.7) + 
  xlab("TCGA cohort") + ylab("AA fraction") + 
  scale_shape_manual(values = c(16, 17)) + 
  theme_classic() + ggtitle("TCGA_motif selection exonic regions") + 
  theme(plot.title = element_text(hjust = 0.5))

total_plot = hans_plot_aa + hans_plot_pval + plot_layout(guides = 'collect')
ggsave("Figures/TCGA_jitter_plot_by_tissue.png", total_plot, width = 10, height = 5)


# plots by signatures
hans_plot_SBS88 = ggplot(meta_TCGA, aes(x = cohort, y =SBS88 , color =  Motif_selection, shape = POLE_somatic)) + 
  geom_jitter(alpha = 0.7) + 
  scale_shape_manual(values = c(16, 17)) +
  xlab("TCGA cohort") + ylab("SBS88 WES contribution") + 
  theme_classic() + ggtitle("TCGA SBS88") + 
  theme(plot.title = element_text(hjust = 0.5))

hans_plot_ID18 =  ggplot(meta_TCGA, aes(x = cohort, y =SBS88 , color =  Motif_selection, shape = POLE_somatic)) + 
  geom_jitter(alpha = 0.7) + 
  scale_shape_manual(values = c(16, 17)) +  xlab("TCGA cohort") + ylab("ID18- contribution") + 
  theme_classic() + ggtitle("ID18 WES contribution") + 
  theme(plot.title = element_text(hjust = 0.5))
sup_plot  = hans_plot_SBS88 + hans_plot_ID18 + plot_layout(guides = 'collect') + 
  plot_annotation(tag_levels = "A")
ggsave("../Manuscript/Figures/Figure_6/Supplementary_TCGA_sigs.pdf", sup_plot, width = 8, height = 5)



p_value_color_plot = ggplot(meta_TCGA, aes(x = AA, y =pval_dinucs , color =  cohort, shape = POLE_somatic)) + 
  geom_point(alpha = 0.6) + 
  geom_hline(yintercept = 3) + geom_vline(xintercept = 0.26) + 
  # scale_color_manual(values = c("lightblue", "#d95f02"), name = "WGS motif selection") + 
  scale_shape_manual(values = c(16, 17)) +  
  xlab("fraction mutations with AA at\n-3-4 position at pks-sites") + ylab("-log10 p-value") + 
  theme_classic() + ggtitle("TCGA_motif selection exonic regions") + 
  theme(plot.title = element_text(hjust = 0.5)) 
p_value_color_plot


TCGA_sig_color = ggplot(meta_TCGA, aes(x = ID18, y =SBS88 , color =  cohort)) + 
  geom_point(alpha = 0.6) + 
  geom_hline(yintercept = 0.05) + geom_vline(xintercept = 0.05) + 
  # scale_color_manual(values = c("lightblue", "#d95f02"), name = "WGS motif selection") + 
  scale_shape_discrete(name = "Pks classification", labels = c("pks negative", "pks+ established", "pks+ new")) + 
  xlab("Id 18 contribution") + ylab("SBS88 contribution") + 
  theme_classic() + ggtitle("TCGA_motif selection exonic regions") + 
  theme(plot.title = element_text(hjust = 0.5)) 
TCGA_color_comparison = TCGA_sig_color + p_value_color_plot + plot_layout(guides = "collect")
ggsave("Figures/TCGA_sig_plot_by_tissue.png", TCGA_color_comparison, width = 10, height = 5)

# insert table HMF WGS versus WES




# Supplementary plot using the barplots by cohort
HMF_h_pval + HMF_h_AA + 
  TCGA_h_pval +   TCGA_h_aa  + 
  
  



# Exploratory analyses below this line :: 
# Plot of exomes bytissue

table(meta_TCGA$cohort)
  p_value_AA_exome = ggplot(meta_TCGA, aes(x = AA, y = pval_dinucs, color  = cohort)) + 
  geom_point(alpha = 0.4) + 
  geom_hline(yintercept = 3) + geom_vline(xintercept = 0.26) + 
  scale_shape_discrete(name = "Pks classification", labels = c("pks negative", "pks+ established", "pks+ new")) + 
  xlab("fraction mutations with AA at\n-3-4 position at pks-sites") + ylab("-log10 p-value") + 
  theme_classic() + scale_color_npg() + 
  ggtitle("TCGA_motif selection exonic regions")
p_value_AA_exome

sigs_exome = ggplot(meta_TCGA, aes(x = ID18, y = SBS88, color = cohort)) + geom_point(alpha = 0.8) + 
  scale_shape_discrete(name = "Pks classification", labels = c("pks negative", "pks+ established", "pks+ new")) + 
  xlab("ID18 relative contribution") + ylab("SBS 88 relative contribution") + 
  theme_classic() + ggtitle("TCGA motif selection exonic regions") + 
  scale_color_npg()
sigs_exome
ggsave("../Manuscript/Figures/Figure_6/TCGA_sig_exome.pdf", sigs_exome, width = 7, height = 4)



# # Plots using different coloring.  - the predictive value of the method is not that great in comparison versus 'real' exomes. 
# # Set filter for samples which are 'likeley POLE mutated - high pval, low fraction of mutations)
# meta_TCGA$high_pval_low_AAfraction = ifelse(meta_TCGA$pval_dinucs > 5 & meta_TCGA$AA < 0.3, 'pole', 'non_pole')
# meta_TCGA$sigpos_10percent_SBS88 = ifelse(meta_TCGA$SBS88 > 0.10 , 'sig', 'non_sig')
# 
# 
# sigs_exome = ggplot(meta_TCGA, aes(x = ID18, y = SBS88, color = high_pval_low_AAfraction)) + geom_point(alpha = 0.4) + 
#   scale_color_manual(values = c("lightblue", "#d95f02"), name = "WGS motif selection") + 
#   scale_shape_discrete(name = "Pks classification", labels = c("pks negative", "pks+ established", "pks+ new")) + 
#   xlab("ID18 relative contribution") + ylab("SBS 88 relative contribution") + 
#   theme_classic() + ggtitle("TCGA motif selection exonic regions")
# sigs_exome
# ggsave("../Manuscript/Figures/Figure_6/TCGA_sig_exome.pdf", sigs_exome, width = 7, height = 4)
# 
# p_value_AA_exome = ggplot(meta_TCGA, aes(x = AA, y = pval_dinucs, color = Motif_selection )) + geom_point(alpha = 0.4) + 
#   geom_hline(yintercept = 3) + geom_vline(xintercept = 0.26) + 
# #  scale_color_manual(values = c("lightblue", "#d95f02"), name = "WGS motif selection") + 
#   scale_shape_discrete(name = "Pks classification", labels = c("pks negative", "pks+ established", "pks+ new")) + 
#   xlab("fraction mutations with AA at\n-3-4 position at pks-sites") + ylab("-log10 p-value") + 
#   theme_classic() + ggtitle("TCGA_motif selection exonic regions")
# p_value_AA_exome
# ggsave("../Manuscript/Figures/Figure_6/TCGA_p_value_AA_exome.pdf", p_value_AA_exome, width = 7, height = 4)

# plot tables of the results using ggtextable
# WGS table 
metadata_noPOLE = metadata[metadata$POLE_somatic == "WT",]
table = table(metadata_noPOLE$rec_tissues, metadata_noPOLE$Motif_selection) %>% as.matrix()
table_motif = cbind(table, round(table[,2]/ rowSums(table[,1:2])*100, 1))
table = table(metadata_noPOLE$rec_tissues, metadata_noPOLE$Sig_classification) %>% as.matrix()
table_sig = cbind(table, round(table[,2]/ rowSums(table[,1:2])*100, 1))

total_table = cbind(table_motif, table_sig)

colnames(total_table)[c(3,6)] = c("Motif class.\npercentage", "Signature class.\npercentage")

table_fig = ggtexttable(total_table, theme = ttheme("light"))
table_fig = table_fig %>% 
  tab_add_title(text = "HMF cohort \ncolibactin-linked fraction WGS", face = "bold", padding = unit(1, "line"))

ggsave("../Manuscript/Figures/Figure_6/HMF_WGS_table.pdf", table_fig, width = 10, height = 5)


# WES table
meta_exome_noPOLE = meta_exome[meta_exome$POLE_somatic == "WT",]
table = table(meta_exome_noPOLE$rec_tissues, meta_exome_noPOLE$MS_WES) %>% as.matrix()
table_motif = cbind(table, round(table[,2]/ rowSums(table[,1:2])*100, 1))

table = table(meta_exome_noPOLE$rec_tissues, meta_exome_noPOLE$WES_sig) %>% as.matrix()
table_sig = cbind(table, round(table[,2]/rowSums(table[,1:2])*100, 1))

total_table = cbind(table_motif, table_sig)
colnames(total_table)[c(3,6)] = c("Motif class.\npercentage", "Signature class.\npercentage")


table_fig = ggtexttable(total_table, theme = ttheme("light"))
table_fig = table_fig %>% 
  tab_add_title(text = "HMF cohort \ncolibactin-linked fraction simulated WES ", face = "bold", padding = unit(1, "line"))

ggsave("../Manuscript/Figures/Figure_6/HMF_WES_table.pdf", table_fig, width = 10, height = 5)

# TCGA table 
meta_TCGA_noPOLE = meta_TCGA[meta_TCGA$POLE_somatic == "WT",]
table = table(meta_TCGA_noPOLE$cohort, meta_TCGA_noPOLE$Motif_selection) %>% as.matrix()
table_motif = cbind(table, round(table[,2]/rowSums(table[,1:2])*100, 1))

table = table(meta_TCGA_noPOLE$cohort, meta_TCGA_noPOLE$Sig_classification) %>% as.matrix()
table_sig = cbind(table, round(table[,2]/ rowSums(table[, 1:2])*100, 1))

total_table = cbind(table_motif, table_sig)
colnames(total_table)[c(3,6)] = c("Motif class.\npercentage", "Signature class.\npercentage")

table_fig = ggtexttable(total_table, theme = ttheme("light"))
table_fig = table_fig %>% 
  tab_add_title(text = "TCGA cohorts \ncolibactin-linked fraction WES ", face = "bold", padding = unit(1, "line"))
ggsave("../Manuscript/Figures/Figure_6/TCGA_WES_table.pdf", table_fig, width = 10, height = 5)



  # get fractions for the TCGA and HMF data
# number of CRC cases in the HMF dataset
sum(metadata$primaryTumorLocation == "Colorectum") # 653
nrow(metadata %>% filter(primaryTumorLocation == "Colorectum" & AA > 0.22 & log_p_wgs > 3)) # 105 WGS
105/656

# Fraction exome
sum(meta_exome$primaryTumorLocation == "Colorectum") # 653
nrow(meta_exome %>% filter(primaryTumorLocation == "Colorectum" & AA_exome > 0.22 & log_p_all_exome > 3)) # 40 Exome

40/653 # 6.1%  of CRC samples using motif selection

sum(meta_TCGA$cohort == "CRC") # 654
nrow(meta_TCGA %>% filter(cohort == "CRC" & AA > 0.22 & pval_dinucs > 3 )) # 654
29/654 # 4.4%  of CRC samples using motif selection



meta_exome$moti



# pROC plots 
library("pROC")


roc_SBS88 = roc(meta_exome$Motif_selection ,meta_exome$SBS88, plot = T)
roc_ID18 = roc(meta_exome$Motif_selection ,meta_exome$ID18, plot = T)
roc_SBS88_WES = roc(meta_exome$Motif_selection ,meta_exome$SBS88_WES)
roc_ID88_WES = roc(meta_exome$Motif_selection ,meta_exome$ID18_WES)
roc_pvalue_WES = roc(meta_exome$Motif_selection, meta_exome$log_p_all_exome)
roc_AA_fraction_WES = roc(meta_exome$Motif_selection, meta_exome$AA_exome)
roc_pvalue = roc(meta_exome$Motif_selection, meta_exome$log_p_wgs)
roc_AA_fraction = roc(meta_exome$Motif_selection, meta_exome$AA)

roc_list = list(roc_SBS88, roc_ID18, roc_SBS88_WES, roc_ID88_WES, roc_pvalue_WES, roc_AA_fraction_WES, 
                roc_AA_fraction, roc_pvalue)

names(roc_list) = c("SBS88 contribution", "ID18 contribution", "SBS88 WES contribution", "ID18 WES contribution", 
                    "pvalue WES contribution", "AA WES fraction", "pvalue WGS", "AA WGS fraction")

library(ggsci)
roc_WES = ggroc(roc_list[1:6]) + scale_color_npg() + theme_classic() +  geom_line(size = 1)
ggsave("../Manuscript/Figures/Figure_6/ROC_WES_sig_values.png", roc_WES, width = 7, height  = 4)
# Add in unfair data 
roc_unfair = ggroc(roc_list) + scale_color_npg() + theme_classic() +  geom_line(size = 1)
ggsave("../Manuscript/Figures/Figure_6//Figure_ROC_all_sig_values.png", roc_unfair, width = 7, height  = 4)


# supplementary plots: indicating the role of POLE 
meta_TCGA$POLE_somatic = factor(meta_TCGA$POLE_somatic, levels = c("WT", "MUT"))
p_value_AA_exome_POLE = ggplot(meta_TCGA, aes(x = AA, y = pval_dinucs, color =  muts, shape = POLE_somatic)) + geom_point(alpha = 0.7) + 
  geom_hline(yintercept = 3) + geom_vline(xintercept = 0.26) + 
    xlab("fraction mutations with AA at\n-3-4 position at pks-sites") + ylab("-log10 p-value") + 
  theme_classic() + ggtitle("TCGA_motif selection exonic regions") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_viridis()
p_value_AA_exome_POLE
ggsave("../Manuscript/Figures/Figure_6/TCGA_p_value_AA_exome.pdf", p_value_AA_exome_ms, width = 7, height = 4)


# Generate compound plot 
compound_plot =
  sigs_exome_HMF + p_value_AA_exome_HMF +  
  p_value_AA_exome_HMF_tissue +   venn + 
  sigs_exome_ms + p_value_AA_exome_ms + 
  roc_WES +
  plot_annotation(tag_levels = "A") + 
  plot_layout(ncol = 2 )


#compound_plot = sigs_exome_ms + p_value_AA_exome_ms + plot_spacer() + plot_spacer() + sigs_exome_HMF + 
# p_value_AA_exome_HMF + venn + plot_annotation(tag_levels = 'A') + plot_layout(nrow = 2, guides = "collect")
ggsave("../Manuscript/Figures/Figure_6/Figure_6v4_R.pdf", compound_plot, width = 12, height = 12)  


# Save supplementary figure X 
supp_X = 
  HMF_h_AA + HMF_h_pval + 
  TCGA_h_aa + TCGA_h_pval

# calculate false positive rate: 
30/(30 + (4822 - 117))

117/

