# Download gnomAD VCF for GRCh37/hg19 (example link, check for the latest version)
wget -c https://example.com/gnomad.genomes.r3.0.sites.vcf.bgz -O gnomad.vcf.bgz
# Index the VCF file
tabix -p vcf gnomad.vcf.bgz

# Download dbSNP VCF (example link, check for the latest version)
wget -c https://example.com/dbSNP.b153.hg19.vcf.gz -O dbSNP.b153.hg19.vcf.gz
# Index the VCF file
tabix -p vcf dbSNP.b153.hg19.vcf.gz

# Download COSMIC VCF (registration and license required)
wget -c https://example.com/cosmic.v90.hg19.vcf.gz -O cosmic.v90.hg19.vcf.gz
# Index the VCF file
tabix -p vcf cosmic.v90.hg19.vcf.gz

# Base Recalibration with gnomAD and dbSNP
gatk BaseRecalibrator \
  -I $OUTPUT_DIR/dedup_reads.bam \
  -R $REF_GENOME \
  --known-sites gnomad.vcf.bgz \
  --known-sites dbSNP.b153.hg19.vcf.gz \
  -O $OUTPUT_DIR/recal_data.table

# Apply BQSR
gatk ApplyBQSR \
  -R $REF_GENOME \
  -I $OUTPUT_DIR/dedup_reads.bam \
  --bqsr-recal-file $OUTPUT_DIR/recal_data.table \
  -O $OUTPUT_DIR/recal_reads.bam

# Somatic Variant Calling with Mutect2 using COSMIC
gatk Mutect2 \
  -R $REF_GENOME \
  -I $OUTPUT_DIR/recal_reads.bam \
  -tumor sample_name \
  --germline-resource gnomad.vcf.bgz \
  --panel-of-normals cosmic.v90.hg19.vcf.gz \
  -O $OUTPUT_DIR/somatic_variants.vcf
  