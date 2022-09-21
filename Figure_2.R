# Nissle Figure 2
# Axel Rosendahl Huber 13-12-2021
library(gtools)
source("E:/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/Scripts/Load_data.R")
setwd("E:/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/")

# select Nissle and Dye I3 exposed clones
I3NIS = c("STE0072I3NISA", "STE0072I3NISB", "STE0072I3NISC")
I3DYE = c("PMC.ORGWT.ISCI3DYEA", "PMC.ORGWT.ISCI3DYEB", "PMC.ORGWT.ISCI3DYEC")

I3NIS_muts = rowSums(mut_mat[, I3NIS])
I3DYE_muts = rowSums((mut_mat[, I3DYE]))

I3_NIS_DYE = data.frame(Nissle = I3NIS_muts, Dye_control = I3DYE_muts, difference = I3NIS_muts - I3DYE_muts)
I3_NIS_DYE[I3_NIS_DYE < 0] = 0


# Boxplots indicating total mutation counts 
dye = colSums(mut_mat[, I3NIS])
nissle = colSums(mut_mat[, I3DYE])
data_frame = data.frame(snv_count = c(dye, nissle), type = rep(c("Nissle", "Control"), each = 3))

comparisons = list( c("Nissle", "Control"))
F2a_sbs_boxplot = ggplot(data_frame, aes(x = type, y = snv_count)) + 
  geom_boxplot() + 
  geom_jitter() +  
  stat_compare_means(comparisons = comparisons, label.y = 650, tip.length = 0.1) +
  ggtitle("singe base substitutions")  +
  theme_BM() + scale_y_continuous(limits = c(0,800))  + 
  theme(plot.title = element_text(hjust = 0.5)) + xlab("")

# indels 
indel_total_counts = c(colSums(indel_counts[,I3NIS]), colSums(indel_counts[,I3DYE])) %>% as.data.frame()
colnames(indel_total_counts) = "indel_count"
indel_total_counts$type = rep(c("Nissle", "Control"), each = 3)
comparisons = list( c("Control", "Nissle"))
F2b_boxplot_indels = ggplot(indel_total_counts, aes(x = type, y = indel_count)) + 
  geom_boxplot() + geom_jitter() +  
  ggtitle("indels") + theme_BM() +
 scale_y_continuous(limits = c(0,80)) + xlab("") +
 stat_compare_means(comparisons = comparisons, label.y = 60, tip.length = 0.2)
  

F2c_sbs_profile = plot_96_profile3(I3_NIS_DYE)

# I3 DYE and I3 Nissle indel analysis
ids_NIS = rowSums(indel_counts[,I3NIS])
ids_DYE = rowSums(indel_counts[,I3DYE])
ids_NIS_DYE = data.frame(Nissle = ids_NIS, Dye_control = ids_DYE, difference = ids_NIS - ids_DYE)
ids_NIS_DYE[ids_NIS_DYE < 0] = 0
F2d_indel_profile = plot_indel_contexts(ids_NIS_DYE) + 
  theme(panel.spacing.x = unit(0, "mm"), legend.position = "none", 
        axis.text.x = element_text(size = 5)) +
  theme_minimal_hgrid() 

# cosine similarity to delta ClbQ & Pks
cossim_sbs = data.frame(Dye_control = I3_NIS_DYE$Dye_control, 
                        Nissle = I3_NIS_DYE$Nissle,
                        Pks = rowSums(mut_mat[, categories == "EWT"])) %>% as.matrix()

F2e_sbs_cosim = cos_sim_matrix(cossim_sbs, cossim_sbs) %>% 
 plot_cosine_heatmap( cluster_cols = F, cluster_rows = F, plot_values = T)

cossim_id = data.frame(Dye_control = ids_NIS_DYE$Dye_control, 
                       Nissle = ids_NIS_DYE$Nissle, 
                       Pks = rowSums(indel_counts[, categories == "EWT"])) %>% as.matrix()

F2f_id_cosim = cos_sim_matrix(cossim_id, cossim_id) %>% 
  plot_cosine_heatmap( cluster_cols = F, plot_values = T, cluster_rows = F)


#  ===== SBS signature re-fitting for each clone ======
sigs_organoid_culture = signatures[, c("SBS1", "SBS5", "SBS18", "SBS88")]
fit_res = fit_to_signatures(mut_mat, as.matrix(sigs_organoid_culture))

fit_res_clones = melt(prop.table(fit_res$contribution, 2))
colnames(fit_res_clones) = c("Signature", "Condition", "Contribution")
fit_res_clones$Exposure = rep(categories, each = 4)
fit_res_clones$Exposure[grepl("I3DYE", fit_res_clones$Condition)] = "Dye_EKO"
fit_res_clones$Exposure[grepl("I3NIS", fit_res_clones$Condition)] = "Nissle"
fit_res_clones = filter(fit_res_clones, grepl("Dye_EKO|EWT|Nissle", Exposure) & Signature == "SBS88")
fit_res_clones$Exposure = factor(fit_res_clones$Exposure, levels = c("Dye_EKO", "Nissle", "EWT"))
levels(fit_res_clones$Exposure) = c("Dye_control", "Nissle", "Pks")

fit_res_clones %>% group_by(Exposure) %>% summarize(max(Contribution))

comparisons = combn(x = unique(fit_res_clones$Exposure %>% as.character()), m = 2) %>% 
  as.data.frame() %>% as.list()

stats = compare_means(Contribution ~ Exposure , fit_res_clones, p.adjust.method = "fdr")
stats$yval = c(8, 0,8, 0.2)

fit_res_clones$Exposure
F2g_sbs_refit = ggplot(fit_res_clones, aes(x = Exposure, y = Contribution, fill = Exposure))  +
  geom_violin(alpha = 0.5) + 
  geom_boxplot(outlier.size = 0, width = 0.15, alpha = 0.6) + 
  geom_jitter(width = 0.15, alpha = 0.8)  + 
  theme_BM() +
  stat_compare_means(comparisons = comparisons[c(3,1,2)]) +
  scale_y_continuous(breaks = seq(0,1.2,0.25), limits =c(-0.005,1.2)) + 
  ylab("Relative SBS88 contribution") + xlab("") +
  theme(legend.position =  "none")
F2g_sbs_refit

# ===== Indel signature refitting ====== 
id_sigs_select = id_signatures[, c("ID1", "ID2", "ID18")]
fit_res_id = fit_to_signatures(indel_counts, as.matrix(id_sigs_select))

fit_res_clones = melt(prop.table(fit_res_id$contribution, 2))
colnames(fit_res_clones) = c("Signature", "Condition", "Contribution")
fit_res_clones$Exposure = rep(categories, each = 3)
fit_res_clones$Exposure[grepl("I3DYE", fit_res_clones$Condition)] = "Dye_EKO"
fit_res_clones$Exposure[grepl("I3NIS", fit_res_clones$Condition)] = "Nissle"
fit_res_clones = filter(fit_res_clones, grepl("Dye_EKO|EWT|Nissle", Exposure) & Signature == "ID18")
fit_res_clones$Exposure = factor(fit_res_clones$Exposure, levels = c("Dye_EKO", "Nissle", "EWT"))
levels(fit_res_clones$Exposure) = c("Dye_control", "Nissle", "Pks")

F2h_id_refit = ggplot(fit_res_clones, aes(x = Exposure, y = Contribution, fill = Exposure))  + 
  geom_violin(alpha = 0.5) + 
  geom_boxplot(outlier.size = 0, width = 0.15, alpha = 0.6) +xlab("") + 
  ylab("Relative ID18 contribution") +
  geom_jitter(width = 0.15, alpha = 0.8)  +
  stat_compare_means(comparisons = comparisons[c(3,1,2)]) +
  theme_BM() + scale_y_continuous(breaks = seq(0,1.2,0.25), limits =c(-0.005,1.2)) + 
  theme(legend.position =  "none")

# patch all plots together
figure_2 = F2a_sbs_boxplot + F2b_boxplot_indels + F2c_sbs_profile + F2d_indel_profile + 
  F2e_sbs_cosim + F2f_id_cosim + 
  F2g_sbs_refit + F2h_id_refit + 
  plot_layout(guides = "collect", byrow = F, nrow = 2,widths = c(0.8, 3.5,1,1.2)) + plot_annotation(tag_levels = "a")
ggsave("../Manuscript/Figures/Figure_2/Figure_2_R.pdf", figure_2, width = 20, height = 9)


