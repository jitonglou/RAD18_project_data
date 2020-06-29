# Codes for generating Figure 2B/2C/2E/2F and Supplementary Table 2

#### Load libraries ####
library(dplyr)
library(magrittr)
library(ggplot2)
library(vcfR)
library(ggpubr) # stat_compare_means
library(scales) # percent
library(gridExtra) # grid.arrange
library(reshape2) # melt


#### Load VCFs ####
getwd()
vcf_dir <- "./RemoveDup/results"
vcf_files <- list.files(vcf_dir, pattern="*_rmdup.vcf")
vcfnames <- c("KO15","KO16","KO18","KO19B","KO19T","KO26",
              "WT1","WT22","WT23","WT6","WT8")

sample_size = data.frame(SAMPLE=vcfnames, 
                         size=c(109705010, # KO15
                                109836164, # KO16
                                110794804, # KO18
                                120628284, # KO19B
                                115026008, # KO19T
                                118077102, # KO26
                                114778038, # WT1
                                108730612, # WT22
                                113104020, # WT23
                                114354364, # WT6
                                109092046 # WT8
                         ))

for (i in 1:length(vcf_files)){
  vcf <- read.table(file.path(vcf_dir,vcf_files[i]), header = F, as.is = c(2:10))
  vcf$SAMPLE <- vcfnames[i]
  if (i==1) {vcfs <- vcf}
  else {vcfs <- rbind(vcfs, vcf)}
}
colnames(vcfs)[1:10] <- c("CHROM","POS","ID","REF","ALT","QUAL",
                          "FILTER","INFO","FORMAT","FORMATVALUE")

#### Count variants by type ####
REFs = c("A","T","A","T","A","T","G","C","G","C","G","C")
ALTs = c("C","G","G","C","T","A","A","T","C","G","T","A")

vcfs_SNV = vcfs %>%
  filter(nchar(REF)==1 & nchar(ALT)==1) %>%
  mutate(type=paste0(REF, ALT)) %>%
  mutate(type_combined=case_when(
    type == "AC" | type == "TG" ~ "ACTG",
    type == "AG" | type == "TC" ~ "AGTC",
    type == "AT" | type == "TA" ~ "ATTA",
    type == "GA" | type == "CT" ~ "GACT",
    type == "GC" | type == "CG" ~ "GCCG",
    type == "GT" | type == "CA" ~ "GTCA"
  )) # count SNVs

vcfs_INDEL = vcfs %>%
  filter(!(nchar(REF)==1 & nchar(ALT)==1)) %>%
  mutate(type = case_when(
    nchar(ALT) > nchar(REF) ~ "INS",
    nchar(ALT) < nchar(REF) ~ "DEL"
  ), LEN = abs(nchar(ALT) - nchar(REF))) %>%
  mutate(type_combined=case_when(
    type == "INS" & LEN <= 4 ~ "INS 4",
    type == "INS" & LEN > 4 ~ "INS 4+",
    type == "DEL" & LEN <= 4 ~ "DEL 4",
    type == "DEL" & LEN > 4 ~ "DEL 4+"
  )) # count INDELs

#### Figure 2B ####
## SNV histogram raw count 
df_hist_SNV = vcfs_SNV %>% 
  filter(!(CHROM == "chr6" & POS > 1e8 & POS < 1.3e8)) %>%
  mutate(Design = gsub("([KOWT]+)[0-9]+[A-z]*","\\1",SAMPLE),
         type_label = as.factor(type_combined))

levels(df_hist_SNV$type_label) = c("A(T)>C(G)","A(T)>G(C)","A(T)>T(A)",
                                   "G(C)>A(T)","G(C)>C(G)","G(C)>T(A)")  

df_hist_SNV %>%
  group_by(type_label, Design, SAMPLE) %>%
  summarise(n=n()) %>%
  group_by(type_label, Design) %>%
  summarise(med=median(n))

## SNV histogram by sample
p_hist_SNV = ggplot(data = df_hist_SNV, mapping = aes(SAMPLE)) +
  geom_bar(aes(fill=Design), width = 0.5) +
  labs(y="Raw number of SNV", x=NULL, fill="") +
  scale_fill_manual(labels = c("Rad18 -/-, 6 samples","Rad18 +/+, 5 samples"),
                    values = c("red", "blue")) +
  theme_bw() +
  theme(legend.position="none", legend.box = "horizontal",
        legend.text = element_text(size = 16, face = "bold")
  ) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 24), # formats of y axis title
        axis.text.y = element_text(size = 24, face = "bold"), # formats of y axis label
        strip.text.x = element_text(size = 22, colour = "black", angle = 0) # formats of facet grid title
  ) + 
  facet_grid(~ type_label)

mypath <- file.path("edit paper","figures","Fig 2","Fig2B_hist_SNV.jpeg")
jpeg(file = mypath, units="in", width=14, height=8.5, res=300) # save high resolution figure
p_hist_SNV
dev.off()

#### Figure 2C ####
## SNV boxplot normalized count
df_box_SNV = df_hist_SNV %>%
  group_by(type_label, SAMPLE) %>%
  summarise(raw_count=n()) %>%
  left_join(sample_size, by="SAMPLE") %>%
  mutate(Design = gsub("([KOWT]+)[0-9]+[A-z]*","\\1",SAMPLE),
         norm_count=raw_count/size*10^8)

p_box_SNV = ggplot(df_box_SNV, aes(x=Design, y=norm_count)) +
  geom_boxplot(aes(fill=Design), width = 0.25) +
  labs(y="Normalized number of SNV", x=NULL, fill=""
  ) +
  scale_fill_manual(labels = c("Rad18 -/-, 6 samples","Rad18 +/+, 5 samples"), 
                    values = c("red", "blue")
  ) +
  theme_bw() +
  theme(legend.position="bottom", legend.box = "horizontal",
        legend.text = element_text(size = 24, face = "bold")
  ) +
  theme(plot.caption = element_text(size = 24),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 24), # formats of y axis title
        axis.text.y = element_text(size = 24, face = "bold"), # formats of y axis label
        strip.text.x = element_text(size = 22, colour = "black", angle = 0) # formats of facet grid title
  ) + 
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("KO", "WT")),
                     label = "p.format", # p.signif, p.format
                     size = 6) + 
  scale_y_log10() +
  facet_grid(~ type_label)

mypath <- file.path("edit paper","figures","Fig 2","Fig2C_box_SNV.jpeg")
jpeg(file = mypath, units="in", width=14, height=8.5, res=300) # save high resolution figure
p_box_SNV
dev.off()

#### Figure 2E ####
## INDEL histogram raw count
df_hist_INDEL = vcfs_INDEL %>% 
  filter(!(CHROM == "chr6" & POS > 1e8 & POS < 1.3e8)) %>%
  mutate(Design = gsub("([KOWT]+)[0-9]+[A-z]*","\\1",SAMPLE),
         type_label = as.factor(type_combined))

levels(df_hist_INDEL$type_label) = c("Deletion (1-4 bp)","Deletion (4+ bp)",
                                     "Insertion (1-4 bp)","Insertion (4+ bp)")  

## INDEL histogram by sample
p_hist_INDEL = ggplot(data = df_hist_INDEL, mapping = aes(SAMPLE)) +
  geom_bar(aes(fill=Design), width = 0.5) +
  labs(y="Raw number of INDEL", x=NULL, fill="") +
  scale_fill_manual(labels = c("Rad18 -/-, 6 samples","Rad18 +/+, 5 samples"), values = c("red", "blue")) +
  theme_bw() +
  theme(legend.position="none", legend.box = "horizontal",
        legend.text = element_text(size = 16, face = "bold")
  ) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 24), # formats of y axis title
        axis.text.y = element_text(size = 24, face = "bold"), # formats of y axis label
        strip.text.x = element_text(size = 22, colour = "black", angle = 0) # formats of facet grid title
  ) + 
  facet_grid(~ type_label)

mypath <- file.path("edit paper","figures","Fig 2","Fig2E_hist_INDEL.jpeg")
jpeg(file = mypath, units="in", width=14, height=8.5, res=300) # save high resolution figure
p_hist_INDEL
dev.off()

#### Figure 2F ####
## INDEL boxplot normalized count
df_box_INDEL = df_hist_INDEL %>%
  group_by(type_label, SAMPLE) %>%
  summarise(raw_count=n()) %>%
  left_join(sample_size, by="SAMPLE") %>%
  mutate(Design = gsub("([KOWT]+)[0-9]+[A-z]*","\\1",SAMPLE),
         norm_count=raw_count/size*10^8)

p_box_INDEL = ggplot(df_box_INDEL, aes(x=Design, y=norm_count)) +
  geom_boxplot(aes(fill=Design), width = 0.25) +
  labs(y="Normalized number of INDEL", x=NULL, fill=""
  ) +
  scale_fill_manual(labels = c("Rad18 -/-, 6 samples","Rad18 +/+, 5 samples"), 
                    values = c("red", "blue")
  ) +
  theme_bw() +
  theme(legend.position="bottom", legend.box = "horizontal",
        legend.text = element_text(size = 24, face = "bold")
  ) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 24), # formats of y axis title
        axis.text.y = element_text(size = 24, face = "bold"), # formats of y axis label
        strip.text.x = element_text(size = 22, colour = "black", angle = 0) # formats of facet grid title
  ) + 
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("KO", "WT")),
                     label = "p.format", # p.signif, p.format
                     size = 6) + 
  scale_y_log10() +
  facet_grid(~ type_label)

mypath <- file.path("edit paper","figures","Fig 2","Fig2F_box_INDEL.jpeg")
jpeg(file = mypath, units="in", width=14, height=8.5, res=300) # save high resolution figure
p_box_INDEL
dev.off()

#### Supplementary Table 2 ####
## Raw SNV counts
SNV_raw = matrix(df_box_SNV$raw_count, ncol = 6, byrow = F)
colnames(SNV_raw) = unique(df_box_SNV$type_label)
rownames(SNV_raw) = unique(df_box_SNV$SAMPLE)

## Raw INDEL counts
INDEL_raw = matrix(df_box_INDEL$raw_count, ncol = 4, byrow = F)
colnames(INDEL_raw) = unique(df_box_INDEL$type_label)
rownames(INDEL_raw) = unique(df_box_INDEL$SAMPLE)

## Normalized SNV counts
SNV_norm = matrix(df_box_SNV$norm_count, ncol = 6, byrow = F)
colnames(SNV_norm) = unique(df_box_SNV$type_label)
rownames(SNV_norm) = unique(df_box_SNV$SAMPLE)

## Normalized INDEL counts
INDEL_norm = matrix(df_box_INDEL$norm_count, ncol = 4, byrow = F)
colnames(INDEL_norm) = unique(df_box_INDEL$type_label)
rownames(INDEL_norm) = unique(df_box_INDEL$SAMPLE)


