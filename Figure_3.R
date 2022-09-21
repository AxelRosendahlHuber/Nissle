# Nissle Figure 3
# Axel Rosendahl Huber 13-12-2021
library(ggseqlogo)
library(ggplot2)
library(gtools)
source("E:/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/Scripts/Load_data.R")
setwd("E:/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/")

# ==== do not use NIS6 data 
# comment section if full table needed (include I6 data)
contexts$Nissle = contexts$Nissle[!grepl("NIS6", names(contexts$Nissle))]

contexts_TN = list()
for (type in names(context_list)) {
  ctx_table= context_list[[type]] %>%
    rbindlist() %>%
    distinct() %>%
    filter(grepl("^T", type))
  contexts_TN[[type]] = ctx_table
}
# =====

####### -3 -4 2bp upstream motif 
ext_context = rbindlist(contexts_TN, idcol = "name") %>% 
  mutate(pos34 =  substr(context, 7,8)) %>% 
  mutate(trinucleotide = factor(trinucleotide, levels = TRIPLETS_96[49:96])) %>% 
  mutate(select =  ifelse(pos34 == "AA", "AA", "other") %>% 
  factor(levels = c("other", "AA")))

# plot profile figure 3e for AA contexts
F3x_AA_context_profile = plot_dinuc_profile(ext_context, dinuc = "AA") + ggtitle(element_blank())

# test p-value for enrichment
list = split(ext_context, ext_context$name)
occurrence_mat = sapply(list, function(x) table(x$select))

fisher.test(occurrence_mat[,c(1,3)], alternative = "greater")$p.value # 2.872848e-210
fisher.test(occurrence_mat[,1:2], alternative = "greater") # 5.402e-06

# perform fisher test on all motif enrichments
list2 = sapply(list, function(x) table(x$pos34))
fisher_table = data.frame(Nissle = rep(NA, 16), Pks = NA, dinucs = dinucs)

for (i in 1:nrow(list2)) {
  select = list2[i,]
  other = colSums(list2[-i,])
  mat = rbind(other, select)
  fisher_table$Nissle[i] = fisher.test(mat[,c(1,2)], alternative = "greater")$p.value
  fisher_table$Pks[i] = fisher.test(mat[,c(1,3)], alternative = "greater")$p.value
}

fisher_table_transform = fisher_table
fisher_table_transform[,1:2] = -log10(fisher_table[,1:2])
fisher_table_transform_m = melt(fisher_table_transform)
F3x_dinc_enrichment = ggplot(fisher_table_transform_m, aes(x  =dinucs, y = value)) + geom_bar(stat = "identity") +
  facet_grid( variable ~. , scales =  "free_y")   + 
  ylab("-log10 pvalue") + theme_classic() + 
  geom_hline(yintercept = -log10(0.05), color = "grey") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7))
F3x_dinc_enrichment

# plot the enrichment of samples at specific motifs: 
l2m = melt(prop.table(list2, 2) , varnames = c("Dinucleotide", "Condition") , value.name = "relative frequency")
l2m$Condition = factor(l2m$Condition, levels = c("Pks", "Nissle", "Negative_control"))

F3_dinuc_frequencies = ggplot(l2m) + geom_bar(aes(x = Dinucleotide, y = `relative frequency` ), stat = "identity") + facet_grid(. ~ Condition , scales = "free_y") + 
  theme_half_open() +  panel_border()  + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6 ), panel.spacing.x = unit(0, "points"))  + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.6), breaks = seq(0, 0.6, 0.2)) + 
  ylab("relative frequency") + xlab("dinucleotide")
F3_dinuc_frequencies

# # plot sequence logo's 
TNctx = list()
for (type in names(contexts_TN)) {
  type_context = contexts_TN[[type]]$context
  TNctx[[type]] = type_context
}

F3a_seqlogo_plots = ggseqlogo(TNctx) + 
  annotate('rect', xmin = 10.5, xmax = 11.5, ymin = 0, ymax = 1.1, alpha = .1, col='black', fill='yellow') +  
  scale_x_continuous(breaks = c(1:21), labels= c(-10:10)) + 
  scale_y_continuous(limits = c(0,1.25), expand = c(0,0)) + 
  xlab("Position") + 
  theme_BM() + 
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1))

#################################################
# Final test: perform pairwise enrichment of nucleotides in each category 
#################################################
# 5 different categories (-4A, -3A, -2 A, -1A, +1T)
contexts = rbindlist(contexts_TN, idcol = "exposure")
contexts_pattern = strsplit(contexts$context,"") %>% unlist() %>% matrix(ncol = 21, byrow = T)
contexts = cbind(contexts_pattern, contexts)

# generate all unique combinations of selected peaks from the Pks motif
enrichments = c("7A", "8A", "9A", "10A", "12T")
nucs_select_2 = combinations(n = 5,r = 2,  enrichments, repeats.allowed = F)  %>% t() %>% as.data.frame() %>%  as.list()
nucs_select_3 = combinations(n = 5, r = 3,v =  enrichments, repeats.allowed = F) %>% t() %>% as.data.frame() %>%  as.list()
nucs_select_4 = combinations(n = 5, r = 4,v =  enrichments, repeats.allowed = F) %>% t() %>% as.data.frame() %>%  as.list()

list_total_enrichments = c(enrichments, nucs_select_2, nucs_select_3, nucs_select_4, list(enrichments))
recode_names = lapply(list_total_enrichments, function(x) dplyr::recode(x, `7A` = "-4A", `8A` = "-3A", `9A` = "-2A", `10A` = "-1A", `12T` = "+1T"))
names(list_total_enrichments) = sapply(recode_names, paste0, collapse = " ")

# combinations of two
temp_table = contexts %>% as.data.frame()
colnames(temp_table)[1:21] = 1:21

test_table = data.frame(names = names(list_total_enrichments),  Pks = rep(NA, length(list_total_enrichments)))

for (i in 1:length(list_total_enrichments)) {
  
  nucs = list_total_enrichments[[i]]
  pos = parse_number(nucs) %>% as.character()
  base = gsub(".*[0-9]", "", nucs)
  
  base_check = temp_table[,pos] == base 
  
  if (length(nucs) > 1 ) { 
    base_check%>% as.matrix()
    idx  = apply(base_check, 1,all)
  } else { idx = base_check}
  
  motif_match = temp_table$exposure[idx] %>% table()
  motif_nomatch = temp_table$exposure[!idx] %>% table()
  
  mat = rbind(motif_match, motif_nomatch)
  
  # fisher exact test for Pks vs control
  test_table[i,2] = fisher.test(mat[2:1,-2], alternative = "greater")$p.value
  
  # fisher exact test for Nissle vs control
  #test_table[i,2]  = fisher.test(mat[2:1,-3], alternative = "greater")$p.value
}

fisher_transf= test_table
fisher_transf$names = factor(fisher_transf$names, levels = fisher_transf$names[order(fisher_transf$Pks, decreasing = F)])
fisher_transf[,2] = -log10(test_table[,2])

legend_mat = matrix('N', nrow = 6, ncol = 31)
colnames(legend_mat) = list_total_enrichments

for (i in 1:ncol(legend_mat)) { 
  name = names(list_total_enrichments)[i]
  index = gsub("[[:alpha:]]|\\+", "", name) %>% strsplit(split = " ") 
  index = as.numeric(index[[1]]) + 5
 # index = index*-1 + 7
  base = gsub("[^[:alpha:]]", "",name) %>% strsplit("")
  legend_mat[index,i] = base[[1]]
}

colnames(legend_mat) = names(list_total_enrichments)
legend_mat[5,] = "T" # set base mutation to T

# reorder matrix
legend_mat_m =  melt(legend_mat, varnames = c("base position", "motif"), value.name = "base")
legend_mat_m$motif = factor(legend_mat_m$motif, levels = levels(fisher_transf$names))
bases = ggplot(legend_mat_m , aes(x = motif, y = `base position`, fill = base)) + geom_tile(color = "white", size = 0.9 ) + theme_void() + 
  #theme(axis.text.y = element_text(size=16)) +
  scale_fill_manual(values = c("green", "grey", 'red')) + 
  theme(legend.position = "none") + 
  scale_y_continuous(breaks = c(1:6))

F3x_position_enrichment = ggplot(fisher_transf, aes(x  =names, y = Pks)) + geom_bar(stat = "identity") + 
  ylab("-log10 pvalue") + theme_classic() + geom_hline(yintercept = -log10(0.05), color = "grey") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0., 220)) + 
  xlab("motif position and bases")
F3x_position_enrichment = F3x_position_enrichment + inset_element(bases, 0, -0.6, 1, 0, ignore_tag = TRUE)
# idea: Make heatmap plot underneath barplot - match the positions and colors 


###### Perform fisher exact test for selected trinucleotide combinations 
# See if further selection of AA motif at selected trinucleotides works better
# do not iterate over all combinations - start wit the most frequently occurring basepair
pks_signature = signatures$SBS88
pks_signature_ordered = signatures$Type[order(pks_signature, decreasing = T)]
pks_signature_ordered = pks_signature_ordered[substr(pks_signature_ordered, 3,3) == "T"] # select only T>N mutations 

# combinations of two
temp_table = contexts %>% as.data.frame()
colnames(temp_table)[1:21] = 1:21
test_table = data.frame(names = pks_signature_ordered,  Pks = rep(NA, length(pks_signature_ordered)))

for (i in 1:length(pks_signature_ordered)) {
  
  trinucs = pks_signature_ordered[1:i]
  trinucs_index = temp_table$trinucleotide %in% trinucs
  AA_index = substr(temp_table$context, 7,8) == "AA"
  idx = trinucs_index & AA_index
  motif_match = temp_table$exposure[idx] %>% table()
  motif_nomatch = temp_table$exposure[!idx] %>% table()
  mat = rbind(motif_match, motif_nomatch)
  
  # fisher exact test for Pks vs control
  test_table[i,2] = fisher.test(mat[2:1,-2], alternative = "greater")$p.value
}

test_table$names = factor(test_table$names, levels = pks_signature_ordered)
test_table[17,] # P-value for EcC vs control-exposed organoids is 6.49193e-237

F3x_stepwise_trinucs = ggplot(test_table, aes(x = names, y = -log10(Pks))) +
  geom_hline(yintercept = max(-log10(test_table$Pks)), linetype = "dashed", color = "darkgrey") + 
  geom_point() + theme_BM() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) + 
  xlab("latest trinucleotide_added")  + ylab("-log10 pvalue")

# conclusion: 17 most frequently occuring trinucleotides is the best value. 

# fisher test for these values for Nissle 
trinucs = pks_signature_ordered[1:17] # select the 19 most commonly mutated samples 
trinucs_index = temp_table$trinucleotide %in% trinucs
AA_index = substr(temp_table$context, 7,8) == "AA"
idx = trinucs_index & AA_index
motif_match = temp_table$exposure[idx] %>% table()
motif_nomatch = temp_table$exposure[!idx] %>% table()
mat = rbind(motif_match, motif_nomatch)
fisher.test(mat[2:1,-3], alternative = "greater") # p-values = 3.116e-07

# step 1: Generate randomly sampled combinations Pks and control data
# sample 10 times for each percentage point. Sample the same number of mutations as Nissle = 983 in total
# step 2: perform fisher exact testtest on presence %>% of pos34 AA presence 
# step 3: perform signature extraction 

set.seed(20210901)

list_fractions = list()
list_pks_sigs = list()
for (fraction_pks in seq(0.0,0.5, 0.05)) {
  
  n_pks = round(983*fraction_pks)
  n_control = round(983*(1-fraction_pks))
  
  data_fraction = data.frame(replicate = 1:10, p_value = rep(NA, 10), odds_ratio = NA, lower_conf = NA, higher_conf = NA)
  
  for (i in 1:10) { 
    
    # pks mutation sampling and counting 
    ctrl_muts = dplyr::slice_sample(contexts_TN$Negative_control, n = n_control, replace = T)
    
    if (fraction_pks > 0 ) {
      pks_muts = dplyr::slice_sample(contexts_TN$Pks, n = n_pks, replace = T)
      total_muts = rbind(pks_muts, ctrl_muts) 
    } else {total_muts = ctrl_muts}
    
    
    test_muts_AA = sum(substr(total_muts$context, 7,8) == "AA" & total_muts$trinucleotide %in% trinucs[1:17])
    test_muts_noAA = nrow(total_muts) - test_muts_AA
    
    
    ctrl_muts_AA = sum(substr(contexts_TN$Negative_control$context, 7,8) == "AA"  & contexts_TN$Negative_control$trinucleotide %in% trinucs[1:17])
    ctrl_muts_noAA = nrow(contexts_TN$Negative_control)- ctrl_muts_AA
    
    # fisher test 
    fmat = matrix(c(test_muts_AA, test_muts_noAA, ctrl_muts_AA, ctrl_muts_noAA), ncol = 2)
    
    ftest = fisher.test(fmat)
    data_fraction[i,] = c(i, ftest$p.value, ftest$estimate, ftest$conf.int[1], ftest$conf.int[2])
    
  }
  
  list_fractions[[as.character(fraction_pks)]] = data_fraction
}

total_fractions = rbindlist(list_fractions, idcol = "fraction")

# calculate  odds ratio for Nissle 
fnis = fisher.test(mat[2:1,-3]) # p-values = 5.219e-05

total_fractions_with_nissle = rbind(total_fractions, data.frame(fraction = "Nissle", replicate = 1, p_value =  fnis$p.value, 
                                                                odds_ratio = fnis$estimate, lower_conf = fnis$conf.int[1], 
                                                                higher_conf = fnis$conf.int[2]))

F3_mixture_plot_nissle =ggplot(total_fractions_with_nissle, aes(x = replicate, y = odds_ratio, ymin = lower_conf, ymax = higher_conf )) +
  geom_pointrange(alpha = 0.6, size = 0.3) + facet_grid( . ~ fraction, scales = "free_x") + 
  theme_half_open() + panel_border(size = 0.5) + geom_hline(yintercept = 1) + geom_hline(yintercept = fnis$estimate) + 
  theme(axis.text.x = element_blank(), panel.spacing.x = unit(0.1, "lines") ) + 
  xlab("fraction pks+ mutations (10 replicates for each condition)") + 
  ylab("odds ratio") +
  theme(strip.text.x = element_text(size = 9), axis.title.x = element_text(size = 12))

# generate preliminary figure 3: 
dinuc_total_muts$AA = dinuc_total_muts$AA + ggtitle(element_blank()) + theme(axis.text.x =  element_text(size = rel(0.8), colour = "black", angle = 90)) + 
  ylab("Number of SBS")

plot_list = F3a_seqlogo_plots +  F3x_position_enrichment + F3x_dinc_enrichment +  F3_dinuc_frequencies + 
  F3x_AA_context_profile +  F3x_stepwise_trinucs + F3_mixture_plot_nissle

layout = "
AAAAAEEEEE
BBBCCEEEEE
DDDFFFGGGG"

total_plot = plot_list + plot_layout(design = layout) +  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = 'bold')) 
total_plot
ggsave("../Manuscript/Figures/Figure_3/Figure_3_R.pdf", total_plot, width = 14, height = 9)
ggsave("../Manuscript/Figures/Figure_3/Figure_3_R.png", total_plot, width = 14, height = 9)



# Figure 3F linear line 



set.seed(20210901)

list_fractions = list()
list_pks_sigs = list()
for (fraction_pks in seq(0.0,0.5, 0.005)) {
  
  n_pks = round(983*fraction_pks)
  n_control = round(983*(1-fraction_pks))
  
  data_fraction = data.frame(replicate = 1:10, p_value = rep(NA, 10), odds_ratio = NA, lower_conf = NA, higher_conf = NA)
  
  for (i in 1:10) { 
    
    # pks mutation sampling and counting 
    ctrl_muts = dplyr::slice_sample(contexts_TN$Negative_control, n = n_control, replace = T)
    
    if (fraction_pks > 0 ) {
      pks_muts = dplyr::slice_sample(contexts_TN$Pks, n = n_pks, replace = T)
      total_muts = rbind(pks_muts, ctrl_muts) 
    } else {total_muts = ctrl_muts}
    
    
    test_muts_AA = sum(substr(total_muts$context, 7,8) == "AA" & total_muts$trinucleotide %in% trinucs[1:17])
    test_muts_noAA = nrow(total_muts) - test_muts_AA
    
    
    ctrl_muts_AA = sum(substr(contexts_TN$Negative_control$context, 7,8) == "AA"  & contexts_TN$Negative_control$trinucleotide %in% trinucs[1:17])
    ctrl_muts_noAA = nrow(contexts_TN$Negative_control)- ctrl_muts_AA
    
    # fisher test 
    fmat = matrix(c(test_muts_AA, test_muts_noAA, ctrl_muts_AA, ctrl_muts_noAA), ncol = 2)
    
    ftest = fisher.test(fmat)
    data_fraction[i,] = c(i, ftest$p.value, ftest$estimate, ftest$conf.int[1], ftest$conf.int[2])
    
  }
  
  list_fractions[[as.character(fraction_pks)]] = data_fraction
}

total_fractions = rbindlist(list_fractions, idcol = "fraction")

# calculate  odds ratio for Nissle 
fnis = fisher.test(mat[2:1,-3]) # p-values = 5.219e-05

total_fractions_with_nissle = rbind(total_fractions, data.frame(fraction = "Nissle", replicate = 1, p_value =  fnis$p.value, 
                                                                odds_ratio = fnis$estimate, lower_conf = fnis$conf.int[1], 
                                                                higher_conf = fnis$conf.int[2]))

F3_mixture_plot_nissle =ggplot(total_fractions_with_nissle, aes(x = fraction, y = odds_ratio, ymin = lower_conf, ymax = higher_conf )) +
  geom_pointrange(alpha = 0.6, size = 0.3)  +
  theme_half_open() + panel_border(size = 0.5) + geom_hline(yintercept = 1) + geom_hline(yintercept = fnis$estimate) + 
  
  set.seed(20210901)

list_fractions = list()
list_pks_sigs = list()
for (fraction_pks in seq(0.0,0.5, 0.01)) {
  
  n_pks = round(983*fraction_pks)
  n_control = round(983*(1-fraction_pks))
  
  data_fraction = data.frame(replicate = 1:10, p_value = rep(NA, 10), odds_ratio = NA, lower_conf = NA, higher_conf = NA)
  
  for (i in 1:1) { 
    
    # pks mutation sampling and counting 
    ctrl_muts = dplyr::slice_sample(contexts_TN$Negative_control, n = n_control, replace = T)
    
    if (fraction_pks > 0 ) {
      pks_muts = dplyr::slice_sample(contexts_TN$Pks, n = n_pks, replace = T)
      total_muts = rbind(pks_muts, ctrl_muts) 
    } else {total_muts = ctrl_muts}
    
    
    test_muts_AA = sum(substr(total_muts$context, 7,8) == "AA" & total_muts$trinucleotide %in% trinucs[1:17])
    test_muts_noAA = nrow(total_muts) - test_muts_AA
    
    
    ctrl_muts_AA = sum(substr(contexts_TN$Negative_control$context, 7,8) == "AA"  & contexts_TN$Negative_control$trinucleotide %in% trinucs[1:17])
    ctrl_muts_noAA = nrow(contexts_TN$Negative_control)- ctrl_muts_AA
    
    # fisher test 
    fmat = matrix(c(test_muts_AA, test_muts_noAA, ctrl_muts_AA, ctrl_muts_noAA), ncol = 2)
    
    ftest = fisher.test(fmat)
    data_fraction[i,] = c(i, ftest$p.value, ftest$estimate, ftest$conf.int[1], ftest$conf.int[2])
    
  }
  
  list_fractions[[as.character(fraction_pks)]] = data_fraction
}

total_fractions = rbindlist(list_fractions, idcol = "fraction")

# calculate  odds ratio for Nissle 
fnis = fisher.test(mat[2:1,-3]) # p-values = 5.219e-05

total_fractions_with_nissle = rbind(total_fractions, data.frame(fraction = "Nissle", replicate = 1, p_value =  fnis$p.value, 
                                                                odds_ratio = fnis$estimate, lower_conf = fnis$conf.int[1], 
                                                                higher_conf = fnis$conf.int[2]))

F3_mixture_plot_nissle =ggplot(total_fractions_with_nissle, aes(x = fraction, y = odds_ratio, ymin = lower_conf, ymax = higher_conf )) +
  geom_pointrange(alpha = 0.6, size = 0.3) +
  theme_half_open() + panel_border(size = 0.5) + geom_hline(yintercept = 1) + geom_hline(yintercept = fnis$estimate) + 
  theme(axis.text.x = element_blank(), panel.spacing.x = unit(0.1, "lines") ) + 
  xlab("fraction pks+ mutations (10 replicates for each condition)") + 
  ylab("odds ratio") +
  theme(strip.text.x = element_text(size = 9), axis.title.x = element_text(size = 12)) +
  xlab("fraction pks+ mutations (10 replicates for each condition)") + 
  ylab("odds ratio") +
  theme(strip.text.x = element_text(size = 9), axis.title.x = element_text(size = 12))
F3_mixture_plot_nissle
#ggsave("../Manuscript/Figures/Figure_3/F3_mixture_plot_range.png", width = 6, height = 6)
ggsave("../Manuscript/Figures/Figure_3/F3_mixture_plot_range_0.01.pdf", width = 6, height = 6)


total_fractions$fraction = as.numeric(total_fractions$fraction)
model = lm(fraction ~ odds_ratio , data = total_fractions) 

estimate_mean = predict(model, newdata = data.frame(odds_ratio = fnis$estimate))
estimate_low = predict(model, newdata = data.frame(odds_ratio = fnis$conf.int[[1]]))
estimate_high = predict(model, newdata = data.frame(odds_ratio = fnis$conf.int[[2]]))
estimate_mean - estimate_low
estimate_high - estimate_mean

 
estimate_low
estimate_high
