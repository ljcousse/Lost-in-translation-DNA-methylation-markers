#!/bin/sh

dir_fastq_1=/share/run231
dir_fastq_2=/share/run239/Melanoma-382642780/FASTQ_Generation_2023-02-27_10_32_07Z-656933458/
dir_work=/data/homes/louisc/Project_KimSmits

dir_genome=$dir_work/Human
dir_cov=$dir_work/cov
dir_fastq=$dir_work/fastq
dir_trimmed=$dir_work/trimmed
dir_mapping=$dir_work/mapping
dir_mapping_bis=$dir_work/mapping_bis
dir_wgbstools=$dir_work/wgbs_tools
dir_output_wgbstools=$dir_work/output_wgbstools

mkdir $dir_genome
mkdir $dir_cov
mkdir $dir_fastq
mkdir $dir_trimmed
mkdir $dir_mapping
mkdir $dir_output_wgbstools

echo $dir_fastq_1
echo $dir_fastq_2
echo $dir_work


########################
# Get reference genome #
########################

cd $dir_genome

# GRCh38
# wget http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
# wget http://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/Homo_sapiens.GRCh38.105.gtf.gz

# GRCh37
wget http://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
wget http://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.toplevel.fa.gz

pigz -d *

#################
# Move raw data #
#################

cd $dir_fastq 

for file in $(find $dir_fastq_1 -mindepth 2)
do
	cp $file .
done

pigz -d *

for file in $(ls)
do 
	mv $file run_231_$file
done

for file in $(find $dir_fastq_2 -mindepth 2)
do
	cp $file .
done

pigz -d *

for file in $(ls | grep -v run_231)
do 
	mv $file run_239_$file
done


############
# Trimming #
############

cd $dir_work

# Trim_galore 
N=20
(
for base_file in $(ls $dir_fastq | cut -d"_" -f 1-4 | uniq)
do 
	((i=i%N)); ((i++==0)) && wait
	trim_galore -q 20 --illumina -o $dir_trimmed --paired --length 15 --retain_unpaired --path_to_cutadapt cutadapt $dir_fastq/$base_file*_R1*.fastq  $dir_fastq/$base_file*_R2*.fastq  &
done
)

##########
# FastQC #
##########

mkdir fastqc
dir_fastqc=$dir_work/fastqc

for file in $(ls $dir_fastq)
do
	fastqc --outdir $dir_fastqc $dir_fastq/$file &
done

mkdir fastqc_trimmed
dir_fastqc_trimmed=$dir_work/fastqc_trimmed

for file in $(ls $dir_trimmed | grep .fq$ | grep -v unpaired)
do
	fastqc --outdir $dir_fastqc_trimmed $dir_trimmed/$file &
done


###########
# Multiqc #
###########

multiqc $dir_fastqc
mv multiqc_report.html multiqc_report_rawdata.html
rm -r multiqc_data

multiqc $dir_fastqc_trimmed
mv multiqc_report.html multiqc_report_trimmeddata.html
rm -r multiqc_data


###########
# Mapping #
###########

cd $dir_work

# Index genome
bismark_genome_preparation $dir_genome

# Mapping Bismark 

N=4
(
for sample in $(ls trimmed/ | grep val | cut -d_ -f 1-4 | uniq)
do 
	((i=i%N)); ((i++==0)) && wait
	bismark -q --non_bs_mm --bowtie2 -p 4 -o $dir_work/mapping --prefix $sample $dir_genome -1 $dir_work/trimmed/$sample\_L001_R1_001_val_1.fq -2 $dir_work/trimmed/$sample\_L001_R2_001_val_2.fq &
done
)

# Stringent mapping
# bismark -q --non_bs_mm --bowtie2 -p 4 --score_min L,0,-0.001 -o $dir_work/mapping --prefix $sample $dir_genome -1 $dir_work/trimmed/$sample\_L001_R1_001_val_1.fq.gz -2 $dir_work/trimmed/$sample\_L001_R2_001_val_2.fq.gz &


###################
# Copy to new map #
###################

cd $dir_work

ls $dir_mapping | grep .bam$ | wc -l

cp -r $dir_mapping $dir_mapping_bis

ls $dir_mapping_bis | grep .bam$ | wc -l

##################
# Sort and Index #
##################

cd $dir_mapping_bis

for file in $(ls | grep ".bam$")
do
	echo $file
	samtools sort $file -o sorted_$file
	rm $file
	samtools index sorted_$file
done

cd $dir_work

ls $dir_mapping_bis | wc -l
ls $dir_mapping_bis | grep .bam$ | wc -l
ls $dir_mapping_bis | grep .bai$ | wc -l

#################################
# Split bam files for amplicons #
#################################

Rscript AmpliconAwareMapping.R

ls $dir_mapping_bis | wc -l
ls $dir_mapping_bis | grep .bam$ | wc -l
ls $dir_mapping_bis | grep .bai$ | wc -l

##############################################
# Resort everything by name & remove indexes #
##############################################

cd $dir_mapping_bis

N=16
(
for file in $(ls | grep ".bai$")
do
	((i=i%N)); ((i++==0)) && wait
	rm $file &
done
) 

ls $dir_mapping_bis | wc -l
ls $dir_mapping_bis | grep .bam$ | wc -l
ls $dir_mapping_bis | grep .bai$ | wc -l

N=16
(
for file in $(ls | grep ".bam$")
do
	((i=i%N)); ((i++==0)) && wait
	samtools sort -n $file -o ${file//_sorted/} &
done
) 

ls $dir_mapping_bis | wc -l
ls $dir_mapping_bis | grep .bam$ | wc -l
ls $dir_mapping_bis | grep .bai$ | wc -l

N=16
(
for file in $(ls | grep ".bam$" | grep sorted)
do
	((i=i%N)); ((i++==0)) && wait
	rm $file &
done
) 

ls $dir_mapping_bis | wc -l
ls $dir_mapping_bis | grep .bam$ | wc -l
ls $dir_mapping_bis | grep .bai$ | wc -l


##########################
# Methylation extraction #
##########################

cd $dir_work

mkdir meth
dir_meth=$dir_work/meth

# Bismark methylation extractor
N=32
(
for file in $(ls mapping_bis/ | grep .bam$)
do
	((i=i%N)); ((i++==0)) && wait
	bismark_methylation_extractor -p --no_overlap --comprehensive --output $dir_meth mapping_bis/$file &
done 
)

sleep 30

ls $dir_meth | wc -l
ls $dir_mapping_bis | grep .bam$ | wc -l
ls $dir_meth | grep CpG | wc -l

cd $dir_meth

# CpG
N=32
(
for file in $(ls | grep CpG)
do 
	((i=i%N)); ((i++==0)) && wait
	#echo $file
	bismark2bedGraph --scaffolds -o $file.bedGraph $file &
done
)

sleep 30

ls $dir_meth | wc -l
ls $dir_mapping_bis | grep .bam$ | wc -l
ls $dir_meth | grep CpG | grep .cov.gz$ | wc -l

N=32
(
for file in $(ls | grep cov.gz$)
do 
	((i=i%N)); ((i++==0)) && wait
	pigz -d $file &
done
)

cd $dir_work

ls $dir_meth | wc -l
ls $dir_meth | grep CpG | grep .cov$ | wc -l
ls $dir_meth | grep CpG | wc -l


ls $dir_meth | grep CpG | grep .cov$ | wc -l
ls $dir_mapping_bis | grep .bam$ | wc -l

###########
# Multiqc #
###########

multiqc $dir_meth $dir_mapping
mv multiqc_report.html multiqc_report_bismark.html
rm -r multiqc_data

##########################################
# Extraction reads over unassigned CpG's #
##########################################


for file in $(ls $dir_mapping | grep .bam$)
do 
	samtools view $dir_mapping/$file | awk '$3==15 && $4>33008581-1 && $4<33011868-1 {print ">"NR"\t"$4"\n"$10}' >> GREM1.fa
done

awk '/^>/{key=$2} {print key, NR, $0}' GREM1.fa | sort -k1,1 -k2,2n | cut -d' ' -f3- > GREM1_sorted.fa

head -2 GREM1_sorted.fa

for file in $(ls $dir_mapping | grep .bam$)
do 
	samtools view $dir_mapping$file | awk '$3==16 && $4>58495163-1 && $4<58499903-1 {print ">"NR"\t"$4"\n"$10}' >> NDRG4.fa
done

awk '/^>/{key=$2} {print key, NR, $0}' NDRG4.fa | sort -k1,1 -k2,2n | cut -d' ' -f3- > NDRG4_sorted.fa

head -2 NDRG4_sorted.fa


#############
# wgbstools #
#############

# initialize reference genome
wgbs_tools/wgbstools init_genome --fasta_path $dir_genome/*.fa hg37_87

# sort bams
dir_mapping_bis_sorted=$dir_work/mapping_bis_sorted 
mkdir $dir_mapping_bis_sorted

N=16
(
for file in $(ls $dir_mapping_bis | grep ".bam$")
do
	((i=i%N)); ((i++==0)) && wait
	samtools sort $dir_mapping_bis/$file -o $dir_mapping_bis_sorted/${file//_sorted/} &
done
) 

# create pat & beta files
wgbs_tools/wgbstools bam2pat --genome hg37_87 $dir_mapping_bis_sorted/*.bam -o $dir_output_wgbstools
wgbs_tools/wgbstools vis segmentation_wgbs/pat/run_231_18_S18.run_231_18_S18_L001_R1_001_val_1_bismark_bt2_pe.pat.gz -r 16:58495134-58495243
