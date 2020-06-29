# Codes for estimating somatic signatures and generating Figure 5A/5B/5C/5D

getwd()

#### Load libraries ####
library(SomaticSignatures)
library(BSgenome.Mmusculus.UCSC.mm10)

#### Data: Load VCFs, preproccessing of the TCGA WES Studies ####
vcfPath <- "D:/projects/mouse genome sequencing/annovar/results_refinedVCF/using_bedtools/results"
vcfFile <- list.files(vcfPath, pattern = '*.vcf$')
vcfdir <- file.path(vcfPath, vcfFile)
vcfnames <- c("KO15","KO16","KO18","KO19B","KO19T","KO26",
              "WT1","WT22","WT23","WT6","WT8")

sample_vr <- vector("list",length(vcfdir))
for (i in 1:length(vcfdir)){
  vr <- readVcfAsVRanges(vcfdir[i],"mm10")
  sampleNames(vr) = vcfnames[i]
  sample_vr[[i]] <- vr[nchar(ref(vr))==1 & nchar(alt(vr))==1 & !(seqnames(vr)=="chr6" & start(ranges(vr))>1e8 & start(ranges(vr))<1.3e8),]
}
dwu_vr <- do.call(c, sample_vr)
head(dwu_vr)

#### Motifs: Extracting the Sequence Context of Somatic Variants ####
dwu_motifs = mutationContext(dwu_vr, BSgenome.Mmusculus.UCSC.mm10) # Only SNV substitutions are currently supported.
head(dwu_motifs)

dwu_mm = motifMatrix(dwu_motifs, normalize = TRUE) # 96*11 in this case
head(round(dwu_mm, 4))

plotMutationSpectrum(dwu_motifs, group = "sampleNames") # observed occurrence of the motifs

#### Assessment: Number of Signatures ####
n_sigs = 2:11 # potential numbers

gof_nmf = assessNumberSignatures(dwu_mm, n_sigs, nReplicates = 5)
# nReplicates:	number of runs should be used for assessing a value of 'nSigs'

plotNumberSignatures(gof_nmf) # base on RSS, choose number of somatic signatures = 4

#### Decomposition: Inferring Somatic Signatures ####
n_sigs = 4 # number of somatic signatures

sigs_nmf = identifySignatures(dwu_mm, n_sigs, nmfDecomposition)


#### Figure 5A ####
library(reshape2)
M = as.data.frame(sigs_nmf@observed) # 96*11
M$motif = rownames(M)

M_long = melt(M, id.vars = "motif", variable.name = "sample")
M_long$context = c("A","C","G","T") %>% rep(each=4) %>%
  paste0(".",(c("A","C","G","T") %>% rep(4))) %>%
  rep(66)

df_M_long = M_long %>%
  separate(col=1, into=c("alteration", "context"), sep=" ", remove=F)

p_my_fig5A = ggplot(df_M_long, aes(x=context, y=value, fill=alteration)) +
  geom_histogram(stat="identity") +
  scale_y_continuous(expand = c(0,0),
                     labels = percent_format(accuracy = 1) # scales lib
  ) +
  scale_fill_manual(values = c("deepskyblue3","black","red3","grey","darkolivegreen3","pink2")) +
  theme_bw() +
  guides(fill = guide_legend(nrow = 1)) +
  theme(legend.position = "top") +
  theme(axis.text.x=element_text(angle=90))+
  labs(x="Motifs", y="Mutation Type Probability") +
  facet_grid(sample~alteration)

mypath <- file.path("edit paper","figures","Fig 5",
                    "Fig5A.jpeg")
jpeg(file = mypath, unit="in", width=13, height=12, res=300) # paper="a4r" and "USr" for rotated ("landscape")
p_my_fig5A
dev.off()



#### Figure 5B ####
W = sigs_nmf@signatures %>% as.data.frame() # 96*4
W$motif = rownames(W)

W_long = melt(W, id.vars = "motif", variable.name = "signature")
W_long$context = c("A","C","G","T") %>% rep(each=4) %>%
  paste0(".",(c("A","C","G","T") %>% rep(4))) %>%
  rep(24)

df_W_long = W_long %>%
  separate(col=1, into=c("alteration", "context"), sep=" ", remove=F)

p_my_fig5B = ggplot(df_W_long, aes(x=context, y=value, fill=alteration)) +
  geom_histogram(stat="identity") +
  scale_y_continuous(expand = c(0,0),
                     labels = percent_format(accuracy = 1) # scales lib
  ) +
  # scale_fill_manual(values = c("blue","skyblue","red","pink")) +
  scale_fill_manual(values = c("deepskyblue3","black","red3","grey","darkolivegreen3","pink2")) +
  theme_bw() +
  guides(fill = guide_legend(nrow = 1)) +
  theme(legend.position = "top") +
  theme(axis.text.x=element_text(angle=90))+
  labs(x="Motifs", y="Mutation Type Contribution") +
  facet_grid(signature~alteration)

mypath <- file.path("edit paper","figures","Fig 5",
                    "Fig5B.jpeg")
jpeg(file = mypath, unit="in", width=13, height=4, res=300) # paper="a4r" and "USr" for rotated ("landscape")
p_my_fig5B
dev.off()

#### Figure 5C ####
H <- sigs_nmf@samples
H_norm = apply(H,1,function(z){
  z/sum(z)
}) %>% t

df_H_norm = as.data.frame(H_norm)
df_H_norm$sample = rownames(df_H_norm)
df_H_norm_long = reshape(df_H_norm, idvar="sample", varying = list(1:4),
                         times = paste0("mS",1:4), timevar = "Signature",
                         v.names = "MRC", direction="long")
p_MRC = ggplot(df_H_norm_long, mapping=aes(sample,Signature,fill=MRC)) +
  geom_tile(colour="black") +
  scale_fill_gradientn(name = "MRC",
                       # guide = guide_colourbar(barwidth = 8, barheight = 1),
                       colours=brewer.pal(9,"YlGnBu"),
                       #colours = terrain.colors(10),
                       limit = c(0,1), breaks = c(0, 0.5, 1), space = "Lab"
  ) +
  theme_minimal()

mypath <- file.path("edit paper","figures","Fig 5",
                    "Fig5C_Mutation_signatures_contribution_heatmap.jpeg")
jpeg(file = mypath, unit="in", width=11, height=11, res=300) # paper="a4r" and "USr" for rotated ("landscape")
p_MRC
dev.off()

#### Figure 5D ####
mypath <- file.path("edit paper","figures","Fig 5",
                    "Fig5D_Mutation_signatures_contribution_histogram.jpeg")
jpeg(file = mypath, unit="in", width=11, height=11, res=300) # paper="a4r" and "USr" for rotated ("landscape")
plotSamples(sigs_nmf) +
  scale_fill_brewer(palette = "Paired")
dev.off()
