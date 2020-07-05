#Load packages
library(MutationalPatterns)
library(tidyverse)
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10)
library(deconstructSigs)
ref_genome <- 'BSgenome.Mmusculus.UCSC.mm10'
##Load vcf files from the folder.
vcf.files <- list.files('~/vcfs',full.names = T)
vcf.files <- vcf.files[which(str_detect(vcf.files,'\\.vcf'))]
samplenames <- str_sub(vcf.files,95,101)
samplenames <- str_remove(samplenames,'\\_.*')
##Read vcfs as GRANGES.
vcfs <- read_vcfs_as_granges(vcf.files, samplenames, 'BSgenome.Mmusculus.UCSC.mm10')
##Extract mutation types.
muts = mutations_from_vcf(vcfs[[1]])
types = mut_type(vcfs[[1]])
##Generate context of the mutations.
context = mut_context(vcfs[[1]], ref_genome)
type_context = type_context(vcfs[[1]], ref_genome)
lapply(type_context, head, 12)
##Summary of the mutation in mouse samples.
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
p1 <- plot_spectrum(type_occurrences)
p2 <- plot_spectrum(type_occurrences, CT = TRUE)
##Plot mutational spectrums.
plot_spectrum(type_occurrences, by = samplenames, CT = TRUE)
##96 mutation types * samples.
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
plot_96_profile(mut_mat,condensed = T)

##Use cosmic 30 cancer signatures (COSMIC2) from package deconstrucSigs.
cancer_signatures <- deconstructSigs::signatures.cosmic
##Calculate cosine similarity of cancer signatures and the mutation signatures of the mouse samples.
cos_sim_samples_signatures = cos_sim_matrix(mut_mat, cancer_signatures)

plot_cosine_heatmap(cos_sim_samples_signatures,
                    col_order = cosmic_order,
                    cluster_rows = TRUE)
##Use NMF method to deconvolute mutation signatures to famous 30 cancer signatures.
fit_res <- fit_to_signatures(mut_mat, cancer_signatures)
##Select the top ten signatures that contribute most.
select <- which(rowSums(fit_res$contribution) > 10)
##Plot relative contribution. 
plot_contribution(fit_res$contribution[select,],
                  cancer_signatures[,select],
                  coord_flip = FALSE,
                  mode = "relative")
##Plot absolute contribution. (As in Fig 5F)
plot_contribution(fit_res$contribution[select,],
                  cancer_signatures[,select],
                  coord_flip = FALSE,
                  mode = "absolute")
##Plot contribution heatmap. (As in Fig 5E)
plot_contribution_heatmap(fit_res$contribution,
                          cluster_samples = TRUE,
                          method = "complete")
