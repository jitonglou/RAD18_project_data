#Load the packages.
library(tidyverse)
library(SomaticSignatures)
library(ggstatsplot)
library(ggpubr)
library(EnvStats)
library(ggbeeswarm)
library(reshape2)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(maftools)
#Data Preparation
#The following data were all downloaded on TCGA portal(https://portal.gdc.cancer.gov/). 
#Merged mRNA data from TCGA Portal - HTSeq-FPKM-UQ
mRNA0
#Keep unique tumor samples.
mRNA <- mRNA0[,-which(duplicated(str_sub(colnames(mRNA0),1,12)))]
mRNA <- mRNA[,which(str_sub(colnames(mRNA0),14,14)=='0')]
#CNV data downloaded from TCGA Portal - Gene Level Copy Number Scores
cnv
#Clinical information downloaded from TCGA Portal - clinical
cli
#MAF file downloaded from TCGA Portal. Use package maftools to read maf file.
maf
##Turn maf file into data frame.
mut <- maf@data
##RAD18 mRNA expression
radmRNA <- as.data.frame(mRNA[which(str_detect(rownames(mRNA),"ENSG00000070950")),])
##Log transformation
radmRNA <- log2(radmRNA+1)
colnames(radmRNA) <- "radnum"
radmRNA$id <- rownames(radmRNA)

##Smoking information
smokerlist <- sup1[,c(which(str_detect(colnames(sup1),"bcr_patient_barcode")),
                      which(str_detect(colnames(sup1),"tobacco_smoking_history_indicator")),
                      which(str_detect(colnames(sup1),"tobacco_smoking_pack_years_smoked")))]
##A pack consists of 20 cigarettes. Turn packs/yr into cigarettes/day.
smokerlist$cigarettes_per_day <- as.numeric(smokerlist$smokerlist$tobacco_smoking_pack_years_smoked)*20/365
##Keep only patients with available precise record.
smokerlist <- smokerlist[which(str_length(smokerlist$tobacco_smoking_history_indicator)==1),]
smokerlist$cigarettes_per_day <- as.numeric(as.character(smokerlist$cigarettes_per_day))
##Current reformed smoker (whose indicators are 3 and 4) were considered as nonsmokers.
e <- smokerlist[which(smokerlist$cigarettes_per_day<1),]
e <- e[-which(e$tobacco_smoking_history_indicator==2),]
##Define the smokers and nonsmokers.
smokerlist$cigar <- NA
smokerlist$cigar[which(smokerlist$tobacco_smoking_history_indicator==1)] <- "Nonsmoker"
smokerlist$cigar[which(smokerlist$bcr_patient_barcode%in%e$bcr_patient_barcode)] <- "Nonsmoker"
smokerlist$cigar[which(smokerlist$cigarettes_per_day>1)] <- "Smoker"
smokerlist$cigar[which(smokerlist$tobacco_smoking_history_indicator==2)] <- "Smoker"
smokerlist$cigarettes_per_day[which(smokerlist$tobacco_smoking_history_indicator==1)] <- NA
##Plot the histogram.
smokerlist$cigarettes_per_day[which(smokerlist$tobacco_smoking_history_indicator==1)] <- 0

##Divide samples by RAD18 high/low.
radmRNA$level <- radmRNA$radnum
radmRNA <- arrange(radmRNA,desc(level))
for (i in 1:nrow(radmRNA)){
  ifelse(radmRNA$level[i]> quantile(radmRNA$radnum,.5),
         radmRNA$level[i] <- "high",
         radmRNA$level[i] <- "low")}
radmRNA <- left_join(radmRNA,smokerlist,by = c("id"="bcr_patient_barcode"))
##Remove samples without smoking information.
radmRNA <- radmRNA[-which(is.na(radmRNA$cigar)),]
radmRNA$rad <- str_c(radmRNA$cigar,"-",radmRNA$level)
##Remove duplicated samples.
radmRNA <- radmRNA[-which(duplicated(radmRNA$id)),]

#WES data process.
mut$Tumor_Sample_Barcode <- str_sub(mut$Tumor_Sample_Barcode,1,12)
mut$Chromosome <- str_replace_all(mut$Chromosome,"chr","")
##SNV study
mut1 <- mut[which(mut$Variant_Type=="SNP"),]
RADEXP <- mut1[which(mut1$Tumor_Sample_Barcode %in% smokerlist$bcr_patient_barcode),]
RADEXP <- inner_join(RADEXP,radmRNA,by = c("Tumor_Sample_Barcode"="id"))
##Reshape mutation data into VRANGES, grouped by RAD18 expression level.
RADEXP <-  VRanges(
  seqnames = RADEXP$Chromosome,
  ranges = IRanges(start = RADEXP$Start_Position, 
                   end = RADEXP$End_Position),
  ref = RADEXP$Reference_Allele,
  alt = RADEXP$Tumor_Seq_Allele2,
  sampleNames = RADEXP$Tumor_Sample_Barcode,
  study = RADEXP$rad)
##Extract the sequence context surrounding SNVs from reference GRCh38.
RADEXP_motifs = mutationContext(RADEXP, BSgenome.Hsapiens.NCBI.GRCh38)

plotMutationSpectrum(RADEXP_motifs,group = "study")

##Reshape motif data frame.
b <- as.data.frame(RADEXP_motifs)
c <- dcast(b, sampleNames+study~alteration)
c <- melt(c)
c$value <- log2(c$value+1)

d <- dcast(b, sampleNames+study~alteration)
d$sum <- rowSums(d[,3:8])
d$sum <- log2(d$sum+1)

##Continuous value of RAD18 expression.
RADEXPc <- mut1[which(mut1$Tumor_Sample_Barcode %in% smokerlist$bcr_patient_barcode),]
RADEXPc <- left_join(RADEXPc,radmRNA,by = c("Tumor_Sample_Barcode"="id"))
RADEXPc <-  VRanges(
  seqnames = RADEXPc$Chromosome,
  ranges = IRanges(start = RADEXPc$Start_Position, 
                   end = RADEXPc$End_Position),
  ref = RADEXPc$Reference_Allele,
  alt = RADEXPc$Tumor_Seq_Allele2,
  sampleNames = RADEXPc$Tumor_Sample_Barcode,
  study = RADEXPc$radnum,
  cigar = RADEXPc$cigar)

RADEXPc_motifs = mutationContext(RADEXPc, BSgenome.Hsapiens.NCBI.GRCh38)
bc <- as.data.frame(RADEXPc_motifs)
bc <- bc[which(!is.na(bc$cigar)),]
cc <- dcast(bc, sampleNames+study+cigar~alteration)
cc <- melt(cc,c("sampleNames","study",'cigar'))
cc$value <- log2(cc$value+1)
dc <- dcast(bc, sampleNames+study+cigar~alteration)
dc$sum <- rowSums(dc[,4:9])
dc$sum <- log2(dc$sum+1)
dc <- dc[-which(is.na(dc$study)),]

#CNV study
radcnv <- cnv[which(str_detect(cnv$Gene.Symbol,"ENSG00000070950")),]
radcnv <- radcnv[,-c(1:3)]
##check if have duplications
table(duplicated(colnames(radcnv)))
radcnv <- data.frame(t(radcnv[1,]),colnames(radcnv))
radcnv$colnames.radcnv. <- str_sub(radcnv$colnames.radcnv.,1,12)
##Annotate 
RADCNV_del <- mut1[which(mut1$Tumor_Sample_Barcode %in% radcnv[which(radcnv$X3314==-1),2]),]
RADCNV_amp <- mut1[which(mut1$Tumor_Sample_Barcode %in% radcnv[which(radcnv$X3314==1),2]),]
RADCNV_wt <- mut1[which(mut1$Tumor_Sample_Barcode %in% radcnv[which(radcnv$X3314==0),2]),]
RADCNV_del$RAD <- "RAD18_DEL"
RADCNV_amp$RAD <- "RAD18_AMP"
RADCNV_wt$RAD <- "RAD18_WT"
RADCNV <- rbind(RADCNV_del,RADCNV_amp,RADCNV_wt)

RADCNV <-  VRanges(
  seqnames = RADCNV$Chromosome,
  ranges = IRanges(start = RADCNV$Start_Position, 
                   end = RADCNV$End_Position),
  ref = RADCNV$Reference_Allele,
  alt = RADCNV$Tumor_Seq_Allele2,
  sampleNames = RADCNV$Tumor_Sample_Barcode,
  study = RADCNV$RAD
)
RADCNV_motifs = mutationContext(RADCNV, BSgenome.Hsapiens.NCBI.GRCh38)
sca_mm <- motifMatrix(RADCNV_motifs,group = "study", normalize = T)
sca_sample <- motifMatrix(RADCNV_motifs, normalize = T)

b2 <- as.data.frame(RADCNV_motifs)
c2 <- dcast(b2, sampleNames+study~alteration)
c2 <- melt(c2)

d2 <- dcast(b2, sampleNames+study~alteration)
d2$sum <- rowSums(d2[,4:9])
d2$sum <- log2(d2$sum+1)
d2 <- d2[-which(is.na(d2$study)),]


#####Plotting figures#####

##Histogram of smoking dose.
ggplot(smokerlist)+
  stat_count(aes(cigarettes_per_day),geom = 'bar',width = 0.3)

##Comparison of RAD18 mRNA expression between tumor tissues and normal tissues
figa <- as.data.frame(mRNA0[which(str_detect(rownames(mRNA0),"ENSG00000070950")),])
colnames(figa) <- 'num'
figa$num <- log2(figa$num+1)
figa$type <- str_sub(colnames(mRNA0),14,14)
figa$type[which(figa$type==1)] <- 'Normal'
figa$type[which(figa$type==0)] <- 'Tumor'
figa$id2 <- colnames(mRNA0)
figa$id <- str_sub(figa$id,1,12)
figa <- arrange(figa,desc(num))
figat <- figa[which(figa$type=='Tumor'),]
which(figa$id2%in%figat$id2[which(duplicated(figat$id))])
figa <- figa[-which(figa$id2%in%figat$id2[which(duplicated(figat$id))]),]

ggplot(figa,aes(type,num))+geom_boxplot(outlier.shape = NA,
                                                width = 0.6)+
  geom_beeswarm(size = 1)+
  stat_compare_means()+
  stat_compare_means(comparisons = list(c('Tumor','Normal')),
                     label = 'p.signif')+
  stat_n_text()+
  labs(y = 'log2(RAD18 mRNA FPKM + 1)')

##Comparison of RAD18 expression among RAD18 CNV status
figb <- radcnv
colnames(figb) <- c('CNV','id')
table(str_sub(radcnv$sample_submitter_id,14,14))
figa$id <- str_sub(colnames(mRNA),1,12)
figa$id <- str_replace_all(figa$id,'\\.','-')
figa <- figa[which(figa$type=='Tumor'),]
figb <- left_join(figb,figa,by = c('id'='id'))
figb$CNV <- as.factor(figb$CNV)
ggplot(figb,aes(CNV,num))+geom_boxplot(outlier.shape = NA, width = 0.6)+
  geom_beeswarm(size = 1)+
  stat_compare_means(comparisons = list(c('0','1'),
                                        c('0','-1'),
                                        c('1','-1')),
                     label = 'p.signif')+
  stat_n_text()+
  labs(y = 'log2(RAD18 mRNA FPKM + 1)')

##Mutation Spectrum Grouping by RAD18 Expression Level and Smoking History
plotMutationSpectrum(RADEXP_motifs,group = "study")

##Mutation Spectrum Grouping by RAD18 CNV status and Smoking History
plotMutationSpectrum(RADCNV_motifs,group = "study")

##Comparison of mutation counts grouped by smoking history and RAD18 expression.
ggplot(c,aes(study,value,fill=study))+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~variable) + #coord_cartesian(ylim = ylim1*1.9)+
  stat_compare_means(comparisons = list(c('Smoker-high','Smoker-low'),
                                        c('Nonsmoker-high','Nonsmoker-low')),
                     label = "p.signif", method='wilcox')+#+
  geom_beeswarm(size = 0.3)+
  stat_n_text()+
  labs(x="Smoking History",y="log2(SNV Counts+1)")+
  scale_fill_manual(values = c("blue",'skyblue','red','pink'))+
  scale_colour_manual(values = c("blue",'skyblue','red','pink'))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

##Comparison of total mutation counts grouped by smoking history and RAD18 expression.
ggplot(d,aes(study,sum,fill=study))+geom_boxplot(outlier.shape = NA)+   
  stat_compare_means(comparisons = 
                       list(c("Smoker-low","Smoker-high"),
                            c("Nonsmoker-high","Nonsmoker-low")),
                     label = "p.signif")+
  labs(x="",y="log2(SNV Counts+1)")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  scale_fill_manual(values = c("blue",'skyblue','red','pink'))+
  stat_n_text()+geom_beeswarm(size = 1)

##Scatter plot of mutation counts and RAD18 mRNA expression.
ggplot(cc,aes(study,value))+geom_jitter()+facet_wrap(~variable+cigar) +
  labs(x="RAD18 mRNA Expression",y="log2(SNV Counts+1)")+
  geom_smooth()+stat_cor()

##Scatter plot of total mutation counts and RAD18 mRNA expression.
ggplot(dc,aes(study,sum))+geom_jitter()+facet_wrap(~cigar) +
  labs(x="RAD18 mRNA Expression",y="log2(SNV Counts+1)")+
  geom_smooth()+stat_cor()

##Comparison of mutation counts grouped by RAD18 CNV status.
ggplot(c2,aes(study,value))+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~variable+cigar) + 
  stat_compare_means(comparisons = list(c("Gain","Neutral"),
                                        c("Gain","Loss"),
                                        c("Neutral","Loss")),
                     label = "p.signif")+
  scale_fill_manual(values = c("skyblue",'red'))+
  geom_beeswarm(size = 1)+
  labs(x="",y="log2(SNV Counts+1)")+stat_n_text()

##Comparison of total mutation counts grouped by RAD18 CNV status.
ggplot(d2,aes(study,sum))+
  facet_wrap(vars(cigar))+
  geom_boxplot(outlier.shape = NA)+
  geom_beeswarm(size = 1)+
  stat_compare_means(comparisons = list(c("Gain","Neutral"),
                                        c("Gain","Loss"),
                                        c("Neutral","Loss")),
                     label = "p.signif")+
  labs(x="",y="log2(SNV Counts+1)")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  stat_n_text()