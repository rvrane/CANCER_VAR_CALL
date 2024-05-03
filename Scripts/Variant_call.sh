#!/bin/bash
#SBATCH --job-name=variant_calling
#SBATCH --output=variant_calling_%j.out
#SBATCH --error=variant_calling_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB

# Load modules or set paths to tools
module load FastQC
module load Trimmomatic
module load BWA
module load Samtools
module load Picard
module load GATK
module load SnpEff
module load ANNOVAR
module load REVEL

# Define file paths
FASTQ1="sample_R1.fastq"
FASTQ2="sample_R2.fastq"
REF_GENOME="reference_genome.fa"
KNOWN_VARIANTS="known_variants.vcf"
OUTPUT_DIR="output"

# Step 1: Quality control
fastqc $FASTQ1 $FASTQ2 -o $OUTPUT_DIR

# Step 2: Trimming
trimmomatic PE $FASTQ1 $FASTQ2 \
  $OUTPUT_DIR/trimmed_R1.fastq $OUTPUT_DIR/trimmed_unpaired_R1.fastq \
  $OUTPUT_DIR/trimmed_R2.fastq $OUTPUT_DIR/trimmed_unpaired_R2.fastq \
  SLIDINGWINDOW:4:20 MINLEN:36

# Step 3: Alignment
bwa mem -t 8 $REF_GENOME $OUTPUT_DIR/trimmed_R1.fastq $OUTPUT_DIR/trimmed_R2.fastq > $OUTPUT_DIR/aligned_reads.sam

# Step 4: Convert SAM to BAM, sort and index
samtools view -bS $OUTPUT_DIR/aligned_reads.sam | samtools sort -o $OUTPUT_DIR/sorted_reads.bam
samtools index $OUTPUT_DIR/sorted_reads.bam

# Step 5: Mark duplicates
java -jar $PICARD MarkDuplicates \
  I=$OUTPUT_DIR/sorted_reads.bam \
  O=$OUTPUT_DIR/dedup_reads.bam \
  M=$OUTPUT_DIR/duplicate_metrics.txt \
  REMOVE_DUPLICATES=true

# Step 6: Base recalibration
gatk BaseRecalibrator \
  -I $OUTPUT_DIR/dedup_reads.bam \
  -R $REF_GENOME \
  --known-sites $KNOWN_VARIANTS \
  -O $OUTPUT_DIR/recal_data.table

gatk ApplyBQSR \
  -R $REF_GENOME \
  -I $OUTPUT_DIR/dedup_reads.bam \
  --bqsr-recal-file $OUTPUT_DIR/recal_data.table \
  -O $OUTPUT_DIR/recal_reads.bam

# Step 7: Variant calling
gatk HaplotypeCaller \
  -R $REF_GENOME \
  -I $OUTPUT_DIR/recal_reads.bam \
  -O $OUTPUT_DIR/raw_variants.vcf

# Step 8: Variant filtering (somatic and germline)
gatk Mutect2 \
  -R $REF_GENOME \
  -I $OUTPUT_DIR/recal_reads.bam \
  -tumor sample_name \
  --germline-resource $KNOWN_VARIANTS \
  -O $OUTPUT_DIR/somatic_variants.vcf

# Step 9: Annotation using SnpEff and ANNOVAR
snpEff ann -v hg19 $OUTPUT_DIR/raw_variants.vcf > $OUTPUT_DIR/annotated_variants.vcf
table_annovar.pl $OUTPUT_DIR/somatic_variants.vcf humandb/ -buildver hg19 -out $OUTPUT_DIR/annovar -remove -protocol refGene,cytoBand,exac03,dbnsfp33a -operation g,r,f,f -nastring . -csvout -polish -xreffile example/gene_xref.txt

# Step 10: Predict variant effects using REVEL
revel --input $OUTPUT_DIR/annotated_variants.vcf --output $OUTPUT_DIR/revel_output.txt

# End of script