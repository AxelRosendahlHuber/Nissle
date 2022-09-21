# Figure 5v2 
# POLE mutation status 

# A: HMF WGS pval vs AA (color-coded by mutation load)
# B: HMF WGS highlight hypermutators by color and point out with a line the ones with POLE mutations (as it is in Axel plot right now)
# C: HMF WGS T>N plots with the light color for all and dark color for AA enrichment
# D: HMF AA enrichment in the top 4 trinucleotides of POLE mut sig
# E: HMF AA enrichment in the rest
library(ggrepel)
library(ggExtra)
library(viridis)
library(ggseqlogo)
source("E:/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Manuscript/Figures/Figure_4/Figure_4_v3.R")
setwd("E:/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/")
source("E:/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/Scripts/Nissle_functions.R")

# get samples with somatic POLE mutations 
POLE_muts = list.files("New_HMF/POLE_hotspot_muts/", full.names = T)
names(POLE_muts) = gsub("_.*", "", basename(POLE_muts))
POLE_changes = sapply(POLE_muts, function(x) strsplit(read.table(x)[,8], "\\|")[[1]][11])
metadata$POLE_somatic = "WT"
metadata$POLE_somatic[match(names(POLE_changes), metadata$sampleId)] = POLE_changes
metadata = metadata %>% mutate(POLE_somatic = replace(POLE_somatic,!grepl("Val411Leu|Pro286Arg|Ser459Phe", POLE_somatic), "WT"))
metadata$log10muts = log10(metadata$n_sbs)


# Figure 5A - mutation numbers HMF cohort by color and POLE mutations
F5A = ggplot(metadata, aes(x = AA, y = log_p_wgs, color = log10muts)) +
  geom_point(alpha = 0.8, size = 1) + 
  geom_hline(yintercept = 3) + geom_vline(xintercept = 0.22) + 
  scale_color_viridis() + 
  geom_text_repel(data = metadata %>% filter(POLE_somatic != "WT"),
                  mapping =  aes(label = POLE_somatic), min.segment.length = .05, color = "black",size = 2) + 
  scale_shape_discrete(name = "Pks classification", labels = c("pks negative", "pks+ established", "pks+ new")) + 
  xlab("fraction mutations with AA at\n-3-4 position at pks-sites") + ylab("-log10 p-value") + theme_classic()
F5A

ggsave("../Manuscript/Figures/Figure_5/Figure_5A.pdf", F5A, width = 4.5, height = 3.5)

# figure 5C: 
# select samples with true pole motifs 
POLE_samples = metadata %>% filter(POLE_somatic != "WT") %>% pull(sampleId)
F5c_96_trinucleotide_matrix = plot_96_profile3(mm_HMF[, POLE_samples])


# get mutation table for exome data # hpc needs to be connected for this plot
POLE_list = vector(mode = "list")
for ( sample in POLE_samples) {
  file = paste0("E:/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/Processed_data/HMF_vcfs/",
                sample, "_10basepair_context_SNVs.txt.gz")
  POLE_list[[sample]]  = vroom::vroom(file)
}

POLE_muts = rbindlist(POLE_list, idcol = "name")
POLE_muts = POLE_muts %>% filter(triplet %in% TRIPLETS_96[49:96]) %>% 
  mutate(pos34 =  substr(context, 7,8)) %>% 
  mutate(trinucleotide = factor(triplet, levels = TRIPLETS_96[49:96])) %>% 
  mutate(select =  ifelse(pos34 == "AA", "AA", "other") %>% 
           factor(levels = c("other", "AA")))
POLE_muts$name = "POLE_Mutated_samples"

F5B = plot_dinuc_profile(context = POLE_muts, dinuc = "AA") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.5)) + 
  scale_alpha_manual(values = c(0.5, 1)) + ggtitle("")
  

# Plot contribution TTT
SBS28_tri = TRIPLETS_96[c(96)]
SBS28_counts = fread("E:/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/New_HMF/SBS_28trinuc_counts.txt")
SBS28_counts = slice(SBS28_counts, match(select_HMF, sample)) # only take the select_HMF samples
SBS28_counts$fraction_AA_pos28 = SBS28_counts$TN_AA_sbs28/SBS28_counts$total_muts
SBS28_counts$fraction_AA_NOTpos28 = SBS28_counts$TN_AA_noSBS28/SBS28_counts$total_muts
SBS_28= cbind(SBS28_counts, metadata)
SBS_28$categories = ifelse(SBS_28$log_p_wgs > 3 & SBS_28$AA > 0.22, "motif+", "motif-")
SBS_28$categories[which(SBS_28$sampleId %in% POLE_samples)] = "POLE_hotspot"

F5D = ggplot(SBS_28, aes(y = fraction_AA_pos28, x = categories, fill = categories)) +
  geom_violin() +
  geom_jitter(width = 0.2, alpha = 0.4, size  = 0.3) +
  theme_classic() + ylab("T[T>G]T positions \n with -3-4 AA motif") +
  theme(legend.position = "none") + xlab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

F5E = ggplot(SBS_28, aes(y = fraction_AA_NOTpos28, x = categories, fill = categories)) +
  geom_violin() +
  geom_jitter(width = 0.2, alpha = 0.4, size = 0.3 ) + 
  theme_classic() + ylab("non T[T>G]T positions \n with -3-4 AA motif") + 
  theme(legend.position = "none") +  xlab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Generate figure 5
figure_5 = (F5A /  F5B)  |  (F5D / F5E ) + 
  plot_layout(widths = c(2,1)) +
  plot_annotation(tag_levels = "A") +
  theme(plot.tag = element_text(face = 'bold'))
ggsave(file = "../Manuscript/Figures/Figure_5/Figure_5_Rv3.pdf",
       figure_5, width = 10, height = 8)

