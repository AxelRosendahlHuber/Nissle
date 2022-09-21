# Make hg38 scripts for to analyze all hg38 microinjection data for pks+ bacteria (CCR E.coli and Nissle E. coli)
library(BSgenome)
library(MutationalPatterns)
library(tidyverse)
library(plyr)
library(data.table)
library(vroom)
library(ggpubr)
library(cowplot)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
source("E:/surfdrive/Shared/pmc_vanboxtel/general/2_Bioinformatics/Scripts/pmc_vanboxtel/Axel_BMseq/Utils.R")
source("E:/surfdrive/Shared/pmc_vanboxtel/general/2_Bioinformatics/Scripts/pmc_vanboxtel/Axel_misc//pks_context_selection_functions.R")
source("E:/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/Scripts/Nissle_functions.R")
setwd("E:/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/")
ref_genome = "BSgenome.Hsapiens.UCSC.hg38"
genes = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)

# --------------------------------
get_context =function(gr, size_context = 10){
  gr = gr[gr$FILTER == "PASS"]
  gr = gr[which(nchar(gr$REF) == 1 )]
  strand = ifelse(gr$REF == "G" | gr$REF == "A", '-', "+")
  start = start(gr)
  ref = as.character(gr$REF)
  alt = unlist(CharacterList(gr$ALT))
  type = paste0(ref, ">", alt)
  chromosome = as.character(seqnames(gr))
  context = getSeq(Hsapiens, chromosome,  start = start - size_context,
                   end =  start + size_context, 
                   strand = strand)
      
  type = mapvalues(type, c("A>C","A>G","A>T","G>C", "G>A","G>T"), 
                   c("T>G", "T>C", "T>A", "C>G", "C>T", "C>A"))
  trinucleotide = paste0(substr(context, 10, 10), "[", type, "]", substr(context, 12,12))
  
  
  context_table = tibble(chr = chromosome, position = start, type = type, strand = strand, 
                         context = as.character(context), trinucleotide = trinucleotide)
  context_table$id = paste0(context_table$chr, "_", context_table$position, "_", context_table$type)
  return(context_table)
}



STE72  = list.files("E:/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Boxtel_General/Data/Mutation_data/SNVs/hg38/Axel_GenoEcoli/STE0072/", recursive = T, full.names = T, pattern = "\\.vcf")
STE76 =  list.files("E:/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Boxtel_General/Data/Mutation_data/SNVs/hg38/Axel_GenoEcoli/STE0076/", recursive = T, full.names = T, pattern = "\\.vcf")

####################################
# Exposed organoids mutation data 
####################################

vcf_files = c(STE72, STE76)
names(vcf_files) = gsub(".*_|.vcf", "", basename(vcf_files))
names(vcf_files) = gsub("-", ".", names(vcf_files))
vcf_files = vcf_files[!grepl("TM0", names(vcf_files))]

vcfs = MutationalPatterns::read_vcfs_as_granges(vcf_files, names(vcf_files), genome =  ref_genome, type = "all")
vcfs_sbs = get_mut_type(vcfs, "snv")
vcfs_indel = get_mut_type(vcfs, "indel")
vcfs_dbs = get_mut_type(vcfs, "dbs") # 
vcfs_mnv =get_mut_type(vcfs, "mbs") # total number of mbs in the dataset is 0

mut_mat = mut_matrix(vcfs_sbs, ref_genome)
mut_mat_s = mut_matrix_stranded(vcfs_sbs, ref_genome, genes)

# get double base substitutions - 
dbs_counts = count_dbs_contexts(vcfs_dbs_exp)
plot_dbs_contexts(rowSums(dbs_counts)) # plot of all double base substitution contexts

# Zero mnvs in all samples
lengths(vcfs_mnv)

# get mutation loads 
indel_loads = lengths(id_contexts) %>% as.data.frame()
colnames(indel_loads) = "total_indels"
indel_loads$in_pks_motif = lengths(id_pks_contexts)
indel_loads$fraction_pksmotif = indel_loads$in_pks_motif/indel_loads$total_indels

# ===== Load different names ===== # 
categories = rep("EWT", ncol(mut_mat))
categories[grepl("EKO", colnames(mut_mat))] = "EKO"
categories[grepl("ETBF", colnames(mut_mat))] = "ETBF"
categories[grepl("NIS", colnames(mut_mat))] = "NIS"
categories[grepl("CDT", colnames(mut_mat))] = "CDT"
categories[grepl("DYE", colnames(mut_mat))] = "DYE"

injections = c(rep(1,3), rep(8,2), rep(5,9), rep(8,1), rep(3,10), 
               rep(8,4), rep(3,3), rep(6,3), rep(3,6))

# Load dinucleotide categories
nucs = c("A", "T", "C", "G")
dinucs = expand.grid(nucs, nucs)
dinucs = paste0(dinucs[,1],dinucs[,2])


####################################
# Extended mutational contexts
####################################
contexts = lapply(vcfs_sbs, get_context)

# save unique context files 
context_list = list(Pks = contexts[grepl("EWT", names(contexts))],
                        Nissle = contexts[grep("I3NIS", names(contexts))],
                        Negative_control = contexts[grepl("EKO|DYE", names(contexts))]) 

for (name in names(context_list)) {
  ctx = context_list[[name]]
  ctx = rbindlist(ctx)
  ctx = dplyr::distinct(ctx)
  ctx = ctx[substr(ctx$type, 1,1) == "T","context"]
  write.table(ctx, file = paste0("Processed_data/Contexts/",name, "_unique_contexts.txt"), quote = F, col.names = F, row.names = F)
}

# indel analysis
id_contexts = MutationalPatterns::get_indel_context(vcfs_indel, ref_genome)
id_pks_contexts = lapply(id_contexts, select_context_indel, type = "Strelka")
indel_counts = count_indel_contexts(id_contexts)

####################################
# Hartwig mutational re-fitting data 
####################################
id_signatures = read_delim("https://cancer.sanger.ac.uk/signatures/documents/451/COSMIC_v3.2_ID_GRCh37.txt") 
signatures = read_delim("https://cancer.sanger.ac.uk/signatures/documents/453/COSMIC_v3.2_SBS_GRCh38.txt") %>% 
  arrange(match(Type, TRIPLETS_96))
artefact_signatures = c("SBS27", "SBS43","SBS45","SBS46", "SBS47","SBS48","SBS49","SBS50","SBS51","SBS52",'SBS53',"SBS54",'SBS55','SBS56','SBS57','SBS58','SBS59','SBS60')
sigs = signatures[,!colnames(signatures) %in% artefact_signatures] # remove signatures marked as artefacts

# load HMF data (exome and whole genome) 
mm_HMF = read.delim("E:/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/New_HMF/20210614_HMF_sbs_matrix_somatics.tsv") 
mm_HMF_exome = fread("E:/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/New_HMF/mm_exome.txt", data.table = F)
select_HMF =  colnames(mm_HMF)[colSums(mm_HMF) > 100 ]  # remove all samples with mutation counts < 100
mm_HMF = mm_HMF[, select_HMF]

ids_HMF = read.delim("E:/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/New_HMF/20210614_HMF_indel_matrix_somatics.tsv")[,select_HMF]
ids_HMF_exome = fread("E:/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/New_HMF/indel_counts_HMF_exome.txt", data.table = F)

# remove all samples with no mutation counts in the exome in both indels and snvs
sbsExSelect = colnames(mm_HMF_exome)[colSums(mm_HMF_exome) > 0]
indelExSelect = colnames(ids_HMF_exome)[colSums(ids_HMF_exome) > 0]
exome_select = intersect(intersect(sbsExSelect, indelExSelect), select_HMF)
mm_HMF_exome = mm_HMF_exome[,exome_select]
ids_HMF_exome = ids_HMF_exome[,exome_select]

# HMF WGS refits 
sbs_fit = fit_to_signatures(mm_HMF %>% as.matrix(), signatures[,-1] %>% as.matrix() ) 
sbs_contri = prop.table(sbs_fit$contribution,2)

id_fit = fit_to_signatures(ids_HMF, id_signatures[,-1] %>% as.matrix())
id_contri = prop.table(id_fit$contribution,2)

# HMF exome refits
sbs_fit_exome = fit_to_signatures(mut_matrix = mm_HMF_exome %>% as.matrix(), signatures = signatures[,-1] %>% as.matrix() ) 
sbs_contri_exome = prop.table(sbs_fit_exome$contribution,2)
sbs_contri_exome %>% colSums() %>% is.na() %>% table()

id_fit_exome = fit_to_signatures(ids_HMF_exome, id_signatures[,-1] %>% as.matrix())
id_contri_exome = prop.table(id_fit_exome$contribution,2)

#########
# Hartwig dinucleotide enrichments
##########
# combine the results with signature extraction information: 
# below a code block so the signature analysis is performed like the analysis in the Nature paper 
metadata = read.delim("E:/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/New_HMF/metadata.tsv")
metadata = metadata[match(select_HMF, metadata$sampleId,),] # reorder metadata rows in the same manner as the dinucleotide matrices
metadata$n_sbs =  colSums(mm_HMF)[match(metadata$sampleId, colnames(mm_HMF))] # get total sbs load
metadata$n_indels =  colSums(ids_HMF)[match(metadata$sampleId, colnames(ids_HMF))] # get total indel load

dinuc_mat = read.delim("E:/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/New_HMF/dinuc_mat_hmf_sbs_peaks.txt")
dinuc_mat = dinuc_mat[select_HMF,]
dinuc_exome = read.delim("E:/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/New_HMF/dinuc_mat_exome_sbs_peaks.txt")
dinuc_exome  = dinuc_exome[exome_select,]
# p_value exome: fisher test against all other samples 
log_p_wgs_dinucs = get_dinuc_enrichment(dinuc_mat)
metadata$log_p_wgs = log_p_wgs_dinucs
metadata$log_p_all_exome =  0 
dinuc_enrichment = get_dinuc_enrichment(dinuc_exome)
rownames(metadata) = metadata$sampleId
metadata[names(dinuc_enrichment), "log_p_all_exome"] = dinuc_enrichment

save.image(file = paste0("Processed_data/Nissle_processed", ".RData"))

