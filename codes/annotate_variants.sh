#### set work directories
sample_dir=/proj/vazirilb/HTSF/141006_UNC17-D00216_0248_AC4P6GANXX
dwu_dir=/proj/cosd_lab/dwuWrite
ref_dir=${dwu_dir}/ref
in_dir=${dwu_dir}/annovar/results_refinedVCF/using_bedtools/vcfs
out_dir=${dwu_dir}/annovar/results_refinedVCF/using_bedtools/results_ens
annovar_dir=${dwu_dir}/annovar
sample_name=RAD18WT8_TTAGGC_L002

#### Path to tools
ANNOVAR_PL=/nas/longleaf/apps/annovar/20180416/annovar
GATK=/nas02/apps/biojars-1.0/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar

#### Download files of mouse genomes (on 11/06/2019)
cd ${annovar_dir}/mousedb
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.p6.genome.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.chr_patch_hapl_scaff.annotation.gff3.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.chr_patch_hapl_scaff.annotation.gtf.gz

gunzip GRCm38.p6.genome.fa.gz
gunzip gencode.vM23.chr_patch_hapl_scaff.annotation.gff3.gz
gunzip gencode.vM23.chr_patch_hapl_scaff.annotation.gtf.gz

#### set execute permissions for gtfToGenePred
chmod a+x gtfToGenePred 

${annovar_dir}/mousedb/gtfToGenePred -genePredExt \
    gencode.vM23.chr_patch_hapl_scaff.annotation.gtf \
	mm10_ensGene_og.txt

#### add a random first field to mm10_ensGene_og.txt using nl (adding the line number to each line)
nl mm10_ensGene_og.txt > mm10_ensGene.txt

#### create database of mouse genomes for ensGene
${ANNOVAR_PL}/retrieve_seq_from_fasta.pl \
    --format ensGene \
	--seqfile GRCm38.p6.genome.fa mm10_ensGene.txt \
	--out mm10_ensGeneMrna.fa

#### annote variants using ANNOVAR
${ANNOVAR_PL}/table_annovar.pl \
    ${in_dir}/${sample_name}_bedtools_chr6adj.vcf --vcfinput \
	${annovar_dir}/mousedb/ \
	--outfile ${out_dir}/${sample_name}_bedtools_final \
	--buildver mm10 --protocol ensGene --operation g

#### Parse the INFO column in the VCF file and extract the information of interested filters
$GATK -T VariantsToTable \
    -R ${ref_dir}/mm10_wgs_genome.fa \
	-V ${out_dir}/${sample_name}_bedtools_final.mm10_multianno.vcf \
	-F CHROM -F POS -F REF -F ALT -F QUAL -F FILTER -F DP \
	-F Func.ensGene -F Gene.ensGene -F GeneDetail.ensGene \
	-F ExonicFunc.ensGene -F AAChange.ensGene \
	-o ${out_dir}/${sample_name}_bedtools_annovar_partresults.table