# Codes for generating Figure 3 and Figure 4

#### Libraries and settings ####
library(biomaRt) # useMart, getBM
library(dplyr)
library(tidyr) # complete
library(magrittr) # %<>%
library(IRanges) # IRanges, findOverlaps, overlapsRanges, reduce
library(ggplot2)
library(ggrepel) # geom_text_repel
library(gridExtra) # grid.arrange
library(limma) # vennCounts, vennDiagram, alias2Symbol
library(ggpubr) # stat_compare_means
library(RColorBrewer)
library(colorspace)
library(ggbeeswarm)

getwd()

#### Import the driver gene list and their alias ####
gene_raw = c("Abcb1a","Alkbh3","Ap1m2","Arid1a","Arid2","Ash1l","Blm","Brpf1","Eed","Ep300",
             "Fat1","Fbxw7","Foxe1","Foxm1","Gli1","Gli2","Gli3","Gltscr1","Gsto2","Hras",
             "Hras1","Kmt2c","Kmt2d","Kras","Kras1","Myh9","Nav3","Nfe2l2","Notch1","Notch2",
             "NRas","Nsd1","Palld","Pik3ca","Ptch1","Ptprm","Rexo1","Rpl36","Rras2","Setd2",
             "Sf3a3","Sfrp2","SHH","Slit3","Smg1","Smo","Sparcl1","Sult1a1","Syne1","Syne2",
             "Tgfbr2","Trp53","Unc13a","Unc13c","Usp9x","Vprbp","Wnk2","Wwc1"
)
load("./edit paper/drafts/driver gene alias.RData") # Gtab
genelist = ifelse(is.na(Gtab[,2]), Gtab[,1], Gtab[,2]) # 58 here, but actually 56 since Kras1 doesn't exist, and Hras1 is equivalent to Hras

#### Import the bed file used to generate VCFs ####
bed = read.table("./RefineRegion/new_regions_sorted.bed",
                 as.is = 2:4)
colnames(bed) = c("Chr","Target_start","Target_end","Description")
levels(bed$Chr) = c(as.character(c(1,10:19,2:9)),"MT","X","Y")

#### Obtain annotations for mouse genes and transcriptions on 2019/11/6 ####
# ensembl 98, mouse M23 (GRCm38.p6)
host="www.ensembl.org"

mmart=useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl",host=host)
feature_list=listAttributes(mmart)
features=c("ensembl_gene_id","external_gene_name","gene_biotype",
           "chromosome_name","start_position","end_position",
           "ensembl_transcript_id","external_transcript_name","transcript_biotype",
           "transcript_start","transcript_end","transcription_start_site",
           "description")
mousedb=getBM(attributes=features,mart=mmart) # database used later

#### Import the annotation results of ANNOVAR ####
anno_dir <- "./annovar/results_refinedVCF/using_bedtools/results_ens"
anno_file <- list.files(anno_dir, pattern="*.csv")

samplenames <- c("KO15","KO16","KO18","KO19B","KO19T","KO26",
                 "WT1","WT22","WT23","WT6","WT8")

for (i in 1:length(anno_file)){
  anno <- read.csv(file.path(anno_dir,anno_file[i]), header = T, as.is = c(2:10))
  anno$Sample <- samplenames[i]
  if (i==1) {annos <- anno}
  else {annos <- rbind(annos, anno)}
}

#### Adjust gene length in the mouse database (mousedb) by the bed file ####
mousedb_gene = unique(mousedb[,c(1:6,13)])
mousedb_gene_adj = data.frame(mousedb_gene, gene_length=0)

for (i in 1:nrow(mousedb_gene_adj)){
  gene_og = IRanges(mousedb_gene_adj$start_position[i], mousedb_gene_adj$end_position[i])
  temp = bed %>%
    filter(Chr==mousedb_gene_adj$chromosome_name[i])
  target = reduce(IRanges(temp$Target_start, temp$Target_end)) # some intervals in the bed file are overlapped, so use reduce to get the union of them
  ov = findOverlaps(gene_og, target) # find overlaps between the i-th gene in the mouse database and intervals in the bed file
  mousedb_gene_adj$gene_length[i] = sum(overlapsRanges(gene_og, target, ov)@width) # sum up the width of these overlapped intervals
}

#### Split the rows with multiple genes in the column "ref.Gene" ####
table(annos$Func.ensGene)
annos_filter = annos %>%
  filter(!(Func.ensGene %in% c("intergenic", "exonic;splicing", "ncRNA_exonic;splicing",
                               "upstream;downstream", "UTR5;UTR3")))
genes_split = strsplit(annos_filter$Gene.ensGene, split=";") # split gene(s) of each row

genes_num = sapply(genes_split, length) # get the number of gene(s) in each row
table(genes_num) # multiple genes in one row

annos_genes_split = annos_filter[rep(which(genes_num>=2), genes_num[which(genes_num>=2)]),] # copy those rows with multiple genes
annos_genes_split$Gene.ensGene = unlist(genes_split[which(genes_num>=2)]) # change gene names
annos_adj = rbind(annos_filter[-which(genes_num>=2),], annos_genes_split) # append adjusted rows and remove original rows

annos_adj$Gene.ensGene_lookup = gsub("*\\.[0-9]+","\\1",annos_adj$Gene.ensGene)
annos_adj = annos_adj %>%
  left_join(mousedb_gene %>% dplyr::select(ensembl_gene_id, external_gene_name),
            by=c("Gene.ensGene_lookup"="ensembl_gene_id")) %>%
  unique()

which(duplicated(annos_adj))
table(annos_adj$ExonicFunc.ensGene)
table(annos_adj$Func.ensGene)


#### SNV: KO/WT Common/Specific genes count ####
colors_fig3 = c("pink","skyblue") 

#### Figure 3A ####
genes_only_type = annos_adj %>%
  filter(Ref!="-" & Alt != "-") %>%
  semi_join(rbind(genes_KOonly,genes_WTonly), by="external_gene_name") %>%
  unique() %>%
  inner_join(df_hist_SNV, by=c("Chr"="CHROM","Start"="POS","Ref"="REF",
                               "Alt"="ALT","Sample"="SAMPLE")) 

genes_only_type_box = genes_only_type %>%
  group_by(type_label, Sample) %>%
  summarise(raw_count=n()) %>%
  left_join(sample_size, by=c("Sample"="SAMPLE")) %>%
  mutate(Design = gsub("([KOWT]+)[0-9]+[A-z]*","\\1",Sample),
         norm_count=raw_count/size*10^8)

pp_genes_only_type_box = ggplot(aes(x=Design, y=norm_count),
                                data = genes_only_type_box %>% filter(type_label=="A(T)>T(A)")) +
  geom_boxplot(aes(fill=Design), width = 0.25) +
  labs(y="Normalized number of INDEL", x=NULL, fill=""
  ) +
  scale_fill_manual(labels = c("Rad18 -/-, 6 samples","Rad18 +/+, 5 samples"), 
                    values = c("red", "blue")
  ) +
  theme_bw() +
  theme(legend.position="none", legend.box = "horizontal",
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
  facet_grid(~ type_label)

mypath <- file.path("edit paper","figures","Fig 3","Fig3A.jpeg")
jpeg(file = mypath, units="in", width=4, height=6.5, res=300) # save high resolution figure
pp_genes_only_type_box +
  scale_fill_manual(values=colors_fig3) +
  geom_beeswarm(size = 2, cex = 3)
dev.off()

#### Figure 3B ####
genes_common_type = annos_adj %>%
  filter(Ref!="-" & Alt != "-") %>%
  semi_join(rbind(genes_KOcommon,genes_WTcommon), by="external_gene_name") %>%
  unique() %>%
  inner_join(df_hist_SNV, by=c("Chr"="CHROM","Start"="POS","Ref"="REF",
                               "Alt"="ALT","Sample"="SAMPLE"))  

genes_common_type_box = genes_common_type %>%
  group_by(type_label, Sample) %>%
  summarise(raw_count=n()) %>%
  left_join(sample_size, by=c("Sample"="SAMPLE")) %>%
  mutate(Design = gsub("([KOWT]+)[0-9]+[A-z]*","\\1",Sample),
         norm_count=raw_count/size*10^8)

pp_genes_common_type_box = ggplot(aes(x=Design, y=norm_count),
                                  data = genes_common_type_box %>% filter(type_label=="A(T)>T(A)")) +
  geom_boxplot(aes(fill=Design), width = 0.25) +
  labs(y="Normalized number of INDEL", x=NULL, fill=""
  ) +
  scale_fill_manual(labels = c("Rad18 -/-, 6 samples","Rad18 +/+, 5 samples"), 
                    values = c("red", "blue")
  ) +
  theme_bw() +
  theme(legend.position="none", legend.box = "horizontal",
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
  facet_grid(~ type_label)

mypath <- file.path("edit paper","figures","Fig 3","Fig3B.jpeg")
jpeg(file = mypath, units="in", width=4, height=6.5, res=300) # save high resolution figure
pp_genes_common_type_box +
  scale_fill_manual(values=colors_fig3) +
  geom_beeswarm(size = 2, cex = 3)
dev.off()



#### INDEL: KO/WT Common/Specific genes count ####
df_hist_INDEL_tmp = df_hist_INDEL # ins and del positions are different
df_hist_INDEL_tmp$POS[df_hist_INDEL_tmp$type=="DEL"] = df_hist_INDEL_tmp$POS[df_hist_INDEL_tmp$type=="DEL"]+1 

#### Figure 3C ####
genes_indel_only_type = annos_adj %>%
  filter(Ref=="-" | Alt == "-") %>%
  semi_join(rbind(genes_indel_KOonly,genes_indel_WTonly), by="external_gene_name") %>%
  unique() %>%
  inner_join(df_hist_INDEL_tmp, by=c("Chr"="CHROM","Start"="POS","Sample"="SAMPLE")) 

genes_indel_only_type_box = genes_indel_only_type %>%
  group_by(type_label, Sample) %>%
  summarise(raw_count=n()) %>%
  left_join(sample_size, by=c("Sample"="SAMPLE")) %>%
  mutate(Design = gsub("([KOWT]+)[0-9]+[A-z]*","\\1",Sample),
         norm_count=raw_count/size*10^8)

pp_genes_indel_only_type_box = ggplot(aes(x=Design, y=norm_count),
                                      # data = genes_indel_only_type_box,
                                      data = genes_indel_only_type_box %>% filter(type_label=="Deletion (4+ bp)")) +
  geom_boxplot(aes(fill=Design), width = 0.25) +
  labs(y="Normalized number of SNV", x=NULL, fill=""
  ) +
  scale_fill_manual(labels = c("Rad18 -/-, 6 samples","Rad18 +/+, 5 samples"), 
                    values = c("red", "blue")
  ) +
  theme_bw() +
  theme(legend.position="none", legend.box = "horizontal",
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
  facet_grid(~ type_label)

mypath <- file.path("edit paper","figures","Fig 3","Fig3C.jpeg")
jpeg(file = mypath, units="in", width=4, height=6.5, res=300) # save high resolution figure
pp_genes_indel_only_type_box +
  scale_fill_manual(values=colors_fig3) +
  geom_beeswarm(size=2)
dev.off()

#### Figure 3D ####
genes_indel_common_type = annos_adj %>%
  filter(Ref=="-" | Alt == "-") %>%
  semi_join(rbind(genes_indel_KOcommon,genes_indel_WTcommon), by="external_gene_name") %>%
  unique() %>%
  inner_join(df_hist_INDEL_tmp, by=c("Chr"="CHROM","Start"="POS","Sample"="SAMPLE")) 

genes_indel_common_type_box = genes_indel_common_type %>%
  group_by(type_label, Sample) %>%
  summarise(raw_count=n()) %>%
  left_join(sample_size, by=c("Sample"="SAMPLE")) %>%
  mutate(Design = gsub("([KOWT]+)[0-9]+[A-z]*","\\1",Sample),
         norm_count=raw_count/size*10^8)

pp_genes_indel_common_type_box = ggplot(aes(x=Design, y=norm_count),
                                        data = genes_indel_common_type_box %>% filter(type_label=="Deletion (4+ bp)")) +
  geom_boxplot(aes(fill=Design), width = 0.25) +
  labs(y="Normalized number of SNV", x=NULL, fill=""
  ) +
  scale_fill_manual(labels = c("Rad18 -/-, 6 samples","Rad18 +/+, 5 samples"), 
                    values = c("red", "blue")
  ) +
  theme_bw() +
  theme(legend.position="none", legend.box = "horizontal",
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
  facet_grid(~ type_label)

mypath <- file.path("edit paper","figures","Fig 3","Fig3D.jpeg")
jpeg(file = mypath, units="in", width=4, height=6.5, res=300) # save high resolution figure
pp_genes_indel_common_type_box +
  # scale_fill_hue(l=80) +
  scale_fill_manual(values=colors_fig3) +
  geom_beeswarm(size=2)
dev.off()

#### Appearance of driver genes in each sample ####
## count of each gene by sample
for (i in 1:length(samplenames)){
  count_tmp = annos_adj %>%
    filter(Sample==samplenames[i] & Ref!="-" & Alt != "-") %>%
    group_by(external_gene_name) %>%
    summarise(n=n())
  colnames(count_tmp)[2] = paste0("n_",samplenames[i])
  if (i==1){gene_count_bysample = count_tmp}
  else {gene_count_bysample = full_join(gene_count_bysample, count_tmp, 
                                        by = "external_gene_name")}
}
gene_count_bysample[is.na(gene_count_bysample)] = 0
gene_count_bysample_adj = gene_count_bysample %>%
  mutate(n_KO = n_KO15+n_KO16+n_KO18+n_KO19B+n_KO19T+n_KO26,
         n_WT = n_WT1+n_WT22+n_WT23+n_WT6+n_WT8) %>%
  mutate(n_total = n_KO + n_WT) %>%
  left_join(mousedb_gene_adj, by = c("external_gene_name"))

## count of each driver gene by sample
driver_count_bysample = annos_adj %>%
  filter(external_gene_name %in% genelist) %>%
  group_by(external_gene_name, Sample) %>%
  summarise(count=n()) %>%
  complete(Sample = samplenames)

#### Driver genes SNV ####
## Classify driver genes (common/specific)
driver_class_bydesign_SNV = gene_count_bysample_adj %>%
  filter(external_gene_name %in% genelist) %>%
  dplyr::select(external_gene_name, n_KO, n_WT) %>%
  arrange(external_gene_name) %>%
  mutate(class = case_when(
    n_KO == 0 ~ "Specific_WT",
    n_WT == 0 ~ "Specific_KO",
    TRUE ~ "Common"
  ))

## Appearance of driver genes annotated with SNV in each sample
driver_count_bysample_SNV = annos_adj %>%
  filter(external_gene_name %in% genelist & Ref!="-" & Alt!="-") %>%
  group_by(external_gene_name, Sample) %>%
  summarise(count=n()) %>%
  complete(Sample = samplenames)
driver_count_bysample_SNV$external_gene_name %>% unique %>% length # number of driver genes annotated with SNV
summary(driver_count_bysample_SNV$count)

df_driver_count_bysample_SNV = driver_count_bysample_SNV %>%
  left_join(driver_class_bydesign_SNV, by="external_gene_name") %>%
  unite("class_name", class, external_gene_name, sep="_", remove=FALSE) %>%
  arrange(class_name)

#### Fig 4B ####
windowsFonts(Times=windowsFont("Times New Roman"))

p_driver_count_bysample_SNV = ggplot(df_driver_count_bysample_SNV, 
                                     aes(x=class_name, y=Sample, fill=count)) +
  geom_tile(color="black") +
  geom_text(aes(x=class_name, y=Sample, label=count), color = "black",
            size = 8/72*25.4) +
  scale_fill_gradient(low = "yellow", high = "red", na.value = "white", 
                      limit = c(1,4), space = "Lab",
                      name="Count") + # scale_*_gradient creates a two colour gradient (low-high), scale_*_gradient2 creates a diverging colour gradient (low-mid-high), scale_*_gradientn creates a n-colour gradient.
  scale_x_discrete(breaks = df_driver_count_bysample_SNV$class_name, 
                   labels = df_driver_count_bysample_SNV$external_gene_name) +
  scale_y_discrete(breaks = df_driver_count_bysample_SNV$Sample) +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 65, vjust = 0.6,
                                   size = 12, hjust = 0.5),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24), # formats of y axis title
        axis.text.y = element_text(size = 12)
  ) + 
  labs(y="Sample", x="Gene name", 
       title = "")

mypath <- file.path("edit paper","figures","Fig 4","Fig4B_driver_gene_appearance_SNV.jpeg")
jpeg(file = mypath, units="in", width=14, height=8.5, res=300) # save high resolution figure
p_driver_count_bysample_SNV
dev.off()

#### Fig 4A ####
driver_count_SNV = annos_adj %>%
  filter(external_gene_name %in% genelist & Ref!="-" & Alt!="-")

df_driver_count_SNV = driver_count_SNV %>%
  left_join(driver_class_bydesign_SNV, by="external_gene_name") %>%
  unite("class_name", class, external_gene_name, sep="_", remove=FALSE) %>%
  arrange(class_name) %>%
  mutate(class_name = factor(class_name, levels = unique(class_name)))

## SNV histogram by variant function and exonic variant function
p_driver_count_SNV_Func_ExFunc_histogram <- ggplot(df_driver_count_SNV, 
                                                   aes(class_name,alpha=ExonicFunc.ensGene)) + 
  geom_bar(aes(fill=Func.ensGene)) + 
  theme_classic() + # classic theme
  scale_x_discrete(breaks = df_driver_count_SNV$class_name, 
                   labels = df_driver_count_SNV$external_gene_name) +
  scale_y_continuous(expand = c(0,0)) +
  scale_alpha_manual(values=c(1, 0.7, 0.4, 0.1), name="Exonic variant function") +
  scale_fill_discrete(name = "Variant function") +
  theme(axis.text.x = element_text(angle = 65, vjust = 0.6,
                                   size = 12, hjust = 0.5),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24), # formats of y axis title
        axis.text.y = element_text(size = 12)
  ) +
  labs(x="", y="Number of mutations", title="") # create the figure

mypath <- file.path("edit paper","figures","Fig 4",
                    "Fig4A_driver_gene_appearance_SNV_FuncExoFunc.jpeg")
jpeg(file = mypath, units="in", width=14, height=8.5, res=300) # save high resolution figure
p_driver_count_SNV_Func_ExFunc_histogram
dev.off()


#### Driver genes INDEL ####
driver_count_bydesign_INDEL = annos_adj %>%
  mutate(Design = gsub("([KOWT])[0-9]+[A-Z]*", "\\1", Sample)) %>%
  filter(external_gene_name %in% genelist) %>%
  filter(Ref=="-" | Alt=="-") %>%
  group_by(external_gene_name, Design) %>%
  summarise(n=n()) %>%
  complete(Design = c("KO","WT"),fill=list(n=0))

driver_class_bydesign_INDEL = driver_count_bydesign_INDEL %>%
  filter(Design == "KO") %>%
  dplyr::select(external_gene_name, n) %>%
  left_join((driver_count_bydesign_INDEL %>%
              filter(Design == "WT") %>%
              dplyr::select(external_gene_name, n)), by="external_gene_name", suffix=c("_KO","_WT")) %>%
  mutate(class = case_when(
    n_KO == 0 ~ "Specific_WT",
    n_WT == 0 ~ "Specific_KO",
    TRUE ~ "Common"
  ))

## Appearance of driver genes annotated with INDEL in each sample
driver_count_bysample_INDEL = annos_adj %>%
  filter(external_gene_name %in% genelist & (Ref=="-" | Alt=="-")) %>%
  group_by(external_gene_name, Sample) %>%
  summarise(count=n()) %>%
  complete(Sample = samplenames)
driver_count_bysample_INDEL$external_gene_name %>% unique %>% length # number of driver genes annotated with INDEL
summary(driver_count_bysample_INDEL$count)

df_driver_count_bysample_INDEL = driver_count_bysample_INDEL %>%
  left_join(driver_class_bydesign_INDEL, by="external_gene_name") %>%
  unite("class_name", class, external_gene_name, sep="_", remove=FALSE) %>%
  arrange(class_name)

#### Fig 4D ####
p_driver_count_bysample_INDEL = ggplot(df_driver_count_bysample_INDEL, 
                                       aes(x=class_name, y=Sample, fill=count)) +
  geom_tile(color="black") +
  geom_text(aes(x=class_name, y=Sample, label=count), color = "black",
            size = 12/72*25.4) +
  scale_fill_gradient(low = "yellow", high = "red", na.value = "white", 
                      limit = c(1,4), space = "Lab",
                      name="Count") + # scale_*_gradient creates a two colour gradient (low-high), scale_*_gradient2 creates a diverging colour gradient (low-mid-high), scale_*_gradientn creates a n-colour gradient.
  scale_x_discrete(breaks = df_driver_count_bysample_INDEL$class_name, 
                   labels = df_driver_count_bysample_INDEL$external_gene_name) +
  scale_y_discrete(breaks = df_driver_count_bysample_INDEL$Sample) +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 65, vjust = 0.6,
                                   size = 12, hjust = 0.5),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24), # formats of y axis title
        axis.text.y = element_text(size = 12)
  ) + 
  labs(y="Sample", x="Gene name", 
       title = "")

mypath <- file.path("edit paper","figures","Fig 4",
                    "Fig4D_driver_gene_appearance_INDEL.jpeg")
jpeg(file = mypath, units="in", width=14, height=8.5, res=300) # save high resolution figure
p_driver_count_bysample_INDEL
dev.off()

#### Fig 4C ####
driver_count_INDEL = annos_adj %>%
  filter(external_gene_name %in% genelist) %>%
  filter(Ref=="-" | Alt=="-")

df_driver_count_INDEL = driver_count_INDEL %>%
  left_join(driver_class_bydesign_INDEL, by="external_gene_name") %>%
  unite("class_name", class, external_gene_name, sep="_", remove=FALSE) %>%
  arrange(class_name) %>%
  mutate(class_name = factor(class_name, levels = unique(class_name)))

## INDEL histogram by variant function and exonic variant function
p_driver_count_INDEL_Func_ExFunc_histogram <- ggplot(df_driver_count_INDEL, 
                                                   aes(class_name,alpha=ExonicFunc.ensGene)) + 
  geom_bar(aes(fill=Func.ensGene)) + 
  theme_classic() + # classic theme
  scale_x_discrete(breaks = df_driver_count_INDEL$class_name, 
                   labels = df_driver_count_INDEL$external_gene_name) +
  scale_y_continuous(expand = c(0,0)) +
  scale_alpha_manual(values=c(1, 0.55, 0.1), name="Exonic variant function") +
  scale_fill_discrete(name = "Variant function") +
  theme(axis.text.x = element_text(angle = 65, vjust = 0.6,
                                   size = 12, hjust = 0.5),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24), # formats of y axis title
        axis.text.y = element_text(size = 12)
  ) +
  labs(x="", y="Number of mutations", title="") # create the figure

mypath <- file.path("edit paper","figures","Fig 4",
                    "Fig4C_driver_gene_appearance_INDEL_FuncExoFunc.jpeg")
jpeg(file = mypath, units="in", width=14, height=8.5, res=300) # save high resolution figure
p_driver_count_INDEL_Func_ExFunc_histogram
dev.off()