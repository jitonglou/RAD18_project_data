#### set work directories
project_dir=/lustre/scr/j/i/jitong
sample_dir=/proj/vazirilb/HTSF/141006_UNC17-D00216_0248_AC4P6GANXX
ref_dir=/proj/cosd_lab/dwuWrite/ref
dwu_dir=/proj/cosd_lab/dwuWrite
sample_name=RAD18WT8_TTAGGC_L002

#### path to tools
BBMAP=/nas02/apps/bbmap-37.00/bbmap/bbmap.sh
SAMTOOLS=/nas02/apps/samtools-1.3.1/bin/samtools
PICARD=/nas02/apps/picard-2.2.4/picard-tools-2.2.4/picard.jar
GATK=/nas02/apps/biojars-1.0/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar
INTERSECTBED=/nas02/apps/bedtools-2.26.0/bin/intersectBed

#### path to external files
REF=${ref_dir}/mm10_wgs_genome.fa
KNOWNSNP=${project_dir}/BQSR/mgp.v3.snps.rsIDdbSNPv137_jitong_addchr_sorted.vcf
KNOWNINDEL=${project_dir}/BQSR/mgp.v3.indels.rsIDdbSNPv137_jitong_addchr_sorted.vcf
BED=${proj_dir}/RefineRegion/new_regions_sorted.bed

#### Convert fastq to sam using BBMap
$BBMAP threads=1 \
  in1=${sample_dir}/${sample_name}_R1_001.fastq.gz \
  in2=${sample_dir}/${sample_name}_R2_001.fastq.gz \
  out=${dwu_dir}/sams_bbmap/${sample_name}_jitong.sam \
  ref=${ref_dir}/mm10_wgs_genome.fa nodisk

#### Convert sam to bam using SAMtools
$SAMTOOLS view -bS \
  ${dwu_dir}/sams_bbmap/${sample_name}_jitong.sam \
  > ${project_dir}/original_bams/${sample_name}_jitong.bam

#### Validate bam using Picard
$PICARD ValidateSamFile \
  I=${project_dir}/original_bams/${sample_name}_jitong.bam \
  MODE=SUMMARY

#### Sort bam using Picard
$PICARD SortSam \
  I=${project_dir}/original_bams/${sample_name}_jitong.bam \
  O=${project_dir}/original_bams/${sample_name}_jitong_sorted.bam \
  SORT_ORDER=coordinate

#### Build bam index using Picard
$PICARD BuildBamIndex \
  I=${project_dir}/original_bams/${sample_name}_jitong_sorted.bam

#### Extract read group information using SAMtools 
flow_cell=`$SAMTOOLS view ${project_dir}/original_bams/${sample_name}_jitong_sorted.bam | head -n 1 | cut -f1 | cut -d ':' -f3`
fc_lane=`$SAMTOOLS view ${project_dir}/original_bams/${sample_name}.new.sorted.bam | head -n 1 | cut -f1 | cut -d ':' -f4`
seq_six=`$SAMTOOLS view ${project_dir}/original_bams/${sample_name}.new.sorted.bam | head -n 1 | cut -f1 | cut -d ':' -f10`

#### Add read group using Picard
$PICARD AddOrReplaceReadGroups \
  I=${project_dir}/original_bams/${sample_name}_jitong_sorted.bam \
  O=${project_dir}/bams/${sample_name}_jitong_sorted_RG.bam \
  RGID=${flow_cell}.${seq_six}.${fc_lane} \
  RGLB=lib_${sample_name} \
  RGPL=Illumina \
  RGPU=${flow_cell}.${seq_six}.${fc_lane} \
  RGSM=${sample_name}

#### Validate bam again using Picard
$PICARD ValidateSamFile \
  I=${project_dir}/bams/${sample_name}_jitong_sorted_RG.bam MODE=SUMMARY


#### Remove duplicates using Picard
$PICARD MarkDuplicates \
  I=${project_dir}/bams/${sample_name}_jitong_sorted_RG.bam \
  O=${project_dir}/bams/${sample_name}_jitong_sorted_RG_rmdup.bam \
  METRICS_FILE=${project_dir}/bams/${sample_name}_jitong_rmdup_metrics.txt \
  CREATE_INDEX=TRUE

#### Indel realignment using GATK
$GATK -T RealignerTargetCreator \
  -R $REF \
  -I ${project_dir}/bams/${sample_name}_jitong_sorted_RG_rmdup.bam \
  -known $KNOWNINDEL \
  -o ${project_dir}/indel_relign/${sample_name}_jitong_realignear.intervals

$GATK -T IndelRealigner \
  -R $REF \
  -I ${project_dir}/bams/${sample_name}_jitong_sorted_RG_rmdup.bam \
  -known $KNOWNINDEL \
  -targetIntervals ${project_dir}/indel_relign/${sample_name}_jitong_realignear.intervals \
  -o ${project_dir}/bams/${sample_name}_jitong_sorted_RG_rmdup_realn.bam


#### Base quality score recalibration using GATK
## Analyze patterns of covariation in the sequence dataset
$GATK -T BaseRecalibrator \
  -R $REF \
  -I ${project_dir}/bams/${sample_name}_jitong_sorted_RG_rmdup_realn.bam \
  -knownSites $KNOWNSNP \
  -knownSites $KNOWNINDEL \
  -o ${project_dir}/BQSR/${sample_name}_jitong_recal_data.table

## Do a second pass to analyze covariation remaining after recalibration
$GATK -T BaseRecalibrator \
  -R $REF \
  -I ${project_dir}/bams/${sample_name}_jitong_sorted_RG_rmdup_realn.bam \
  -knownSites $KNOWNSNP \
  -knownSites $KNOWNINDEL \
  -BQSR ${project_dir}/BQSR/${sample_name}_jitong_recal_data.table \
  -o ${project_dir}/BQSR/${sample_name}_jitong_post_recal_data.table 

## Generate before/after plots
$GATK -T AnalyzeCovariates \
  -R $REF \
  -before ${project_dir}/BQSR/${sample_name}_jitong_recal_data.table \
  -after ${project_dir}/BQSR/${sample_name}_jitong_post_recal_data.table \
  -plots ${project_dir}/BQSR/${sample_name}_jitong_recalibration_plots.pdf

## Apply the recalibration to your sequence data
$GATK -T PrintReads \
  -R $REF \
  -I ${project_dir}/bams/${sample_name}_jitong_sorted_RG_rmdup_realn.bam \
  -BQSR ${project_dir}/BQSR/${sample_name}_jitong_recal_data.table \
  -o ${project_dir}/bams/${sample_name}_jitong_sorted_RG_rmdup_realn_recal.bam


#### Variant calling using GATK haplotype caller
$GATK -T HaplotypeCaller \
  -R $REF \
  -I ${project_dir}/bams/${sample_name}_jitong_sorted_RG_rmdup_realn_recal.bam \
  -o ${project_dir}/HCaller/${sample_name}_jitong_HCall.vcf

##### Targeted variant calling using GATK haplotype caller and bed file
$INTERSECTBED \
  -a ${project_dir}/HCaller/${sample_name}_jitong_HCall.vcf \
  -b $BED \
  > ${proj_dir}/HCaller/RefinedVCF/vcf_bedtools/${sample_name}_jitong_HCall_bedtools.vcf
  
