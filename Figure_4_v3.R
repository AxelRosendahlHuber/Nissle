# Figure 4 results
library(ggsci)
setwd("E:/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/")
source("E:/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/Scripts/Load_data.R")
# get relative AA values for each sample (dinuc mat has the same order)
metadata$AA = prop.table(as.matrix(dinuc_mat),1)[,"AA"]
metadata$AA_exome[match(rownames(dinuc_exome), metadata$sampleId)] = prop.table(as.matrix(dinuc_exome),1)[,"AA"]

# fisher test against all other samples 
dinuc_mat_order = dinuc_mat[metadata$sampleId,]

# SBS re-fitting
metadata$SBS88 = sbs_contri["SBS88", match(metadata$sampleId, colnames(sbs_contri))] %>% as.numeric()
metadata$MMR = colSums(sbs_contri[c("SBS6", "SBS14", "SBS15", "SBS20", "SBS21", "SBS26", "SBS44"),  match(metadata$sampleId, colnames(sbs_contri))])
metadata$ID18 = id_contri[ "ID18", match(metadata$sampleId, colnames(id_contri))] %>% as.numeric()
metadata$Sig_classification = ifelse(metadata$SBS88 > 0.05 & metadata$ID18 > 0.05, "sig+", "sig-")
metadata$Motif_selection = ifelse(metadata$log_p_wgs > 3 & metadata$AA > 0.22, "motif+","motif-")
metadata$Motif_Sig = paste0(metadata$Motif_selection, "/", metadata$Sig_classification)
table(metadata$Motif_Sig)
table(metadata$primaryTumorLocation)

# select tissues with > 1 case of colibactin presence 
tissues = table(metadata$Motif_selection, metadata$primaryTumorLocation)
recurrent_tissues = colnames(tissues)[tissues[2,] >= 1]
metadata$rec_tissues = metadata$primaryTumorLocation
metadata$rec_tissues[!metadata$primaryTumorLocation %in% recurrent_tissues] = "Other"
table(metadata$rec_tissues)
table(metadata$primaryTumorLocation[metadata$Motif_selection == "motif+"])
table(metadata$primaryTumorLocation[metadata$Motif_Sig == "motif-/sig+"])

metadata[metadata$Motif_Sig == "motif-/sig+", ]


# p-value plot: 
p_value_AA = ggplot(metadata, aes(x = AA, y = log_p_wgs, color = rec_tissues)) + geom_point(alpha = 0.7) + 
  geom_hline(yintercept = 3) + geom_vline(xintercept = 0.22) + 
  scale_color_manual(values = c("lightblue", "#d95f02", "#1b9e77",  "#bdbdbd", "#7570b3", "black", "maroon", 'pink', 'violet'), name = "Primary cancer origin") + 
  scale_shape_discrete(name = "Pks classification", labels = c("pks negative", "pks+ established", "pks+ new")) + 
  xlab("fraction mutations with AA at\n-3-4 position at pks-sites") + ylab("-log10 p-value") + 
  theme_classic() 

# p-value plot: 
p_value_AA_tissue_co = ggplot(metadata, aes(x = AA, y = log_p_wgs, color = Motif_Sig)) + geom_point(alpha = 0.7) + 
  geom_hline(yintercept = 3) + geom_vline(xintercept = 0.22) + 
  scale_shape_discrete(name = "Pks classification", labels = c("pks negative", "pks+ established", "pks+ new")) + 
  xlab("fraction mutations with AA at\n-3-4 position at pks-sites") + ylab("-log10 p-value") + 
  theme_classic() 
ggsave("../Manuscript/Figures/Figure_4/p_value_AA_tissue_co.pdf", p_value_AA_tissue_co, width = 5.5, height = 4.5)


# confusion matrix 
metadata$pks_fit = "pks-"
metadata$pks_fit[metadata$Motif_selection == "motif+" & metadata$Sig_classification == "sig+"] = "colibactin_established"
metadata$pks_fit[metadata$Motif_selection == "motif+" & metadata$Sig_classification == "sig-"] = "colibactin_new"
conf_mat = table(metadata$Sig_classification, metadata$Motif_selection)

rownames(conf_mat) = c("Sig-","Sig+"); colnames(conf_mat) =  c("Motif-", "Motif+") # set row and column names

# get the cancer types which are motif+ 
metadata %>% filter(Motif_selection == "motif+") %>% pull(primaryTumorLocation) %>% table()
metadata %>% filter(Sig_classification == "sig+") %>% pull(primaryTumorLocation) %>% table()
metadata %>% filter(Motif_Sig == "motif-/sig+") %>% pull(primaryTumorLocation) %>% table()


conf_plot = ggtexttable(conf_mat, theme = ttheme("blank")) %>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)


overview_table = table1::table1( ~ Motif_selection + Sig_classification | rec_tissues , data = metadata)
png("../Manuscript/Figures/Figure_4/overvieuw_table.png")
overview_table
dev.off()

# Todo: Add in 4 colors for the  four different categories: Motif only, Sig only, Motif and sig, and none. 
# Keep coloring across the different plots the same.  

library(table1)
library(grid)
library(gtable)
library(gridExtra)

# correlation plots 
# 4 count classification:



alph = 0.5
SBS_AA_plot = ggplot(metadata , aes(x = AA, y = SBS88, color = Motif_Sig)) + 
  geom_point(alpha = alph) + theme_classic() + 
  theme(legend.position = "none") + xlab("fraction mutations with AA at\n-3-4 position at pks-sites") + 
  ylab("SBS88 contribution")

ID_AA_plot = ggplot(metadata , aes(x = AA, y = ID18, color = Motif_Sig)) + 
  geom_point(alpha = alph) +  theme_classic() +
  theme(legend.position = "none") + xlab("fraction mutations with AA at\n-3-4 position at pks-sites") + 
  ylab("ID18 contribution")

SBS_ID_plot = ggplot(filter(metadata, rec_tissues != "Colorectum") ,
                     aes(x = ID18, y = SBS88, color = Motif_Sig)) + 
  geom_point(alpha = alph) + scale_shape_discrete(name = "Primary cancer origin")  + 
  geom_point(data = filter(metadata, rec_tissues == "Colorectum") ,
             aes(x = ID18, y = SBS88, color = Motif_Sig),
             alpha = alph) + 
  theme_classic() + ylab("SBS88 contribution") + xlab("ID18 contribution") + 
    scale_color_discrete(name = "Pks classification")

# generate figure 4 and supplementary figure
fig4_upper = p_value_AA +  conf_plot + patchwork::plot_spacer()
fig4_lower = SBS_ID_plot + SBS_AA_plot + ID_AA_plot 

fig4 = fig4_upper / fig4_lower 

ggsave("../Manuscript/Figures/Figure_4/Figure_4v3_R.pdf", fig4, width = 10, height = 7)
ggsave("../Manuscript/Figures/Figure_4/Figure_4v3_R.png", fig4, width = 10, height = 7)

