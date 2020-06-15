#!/bin/bash

# get chr1 and simplify
for i in {21..40}
do
cd /lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs
bcftools annotate -x INFO,^FORMAT/GT,FORMAT/PS,FORMAT/DP -r chr1 5739-JL-00"$i".vcf.gz -o 5739-JL-00"$i"_chr1.vcf.gz -O z &
done

# tabix
for i in {21..40}
do
cd /lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs
tabix -p vcf "5739-JL-00"$i"_chr1.vcf.gz" &
done

# run after previous has been finished
vcf-merge \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0021_chr1.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0022_chr1.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0023_chr1.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0024_chr1.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0025_chr1.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0026_chr1.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0027_chr1.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0028_chr1.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0029_chr1.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0030_chr1.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0031_chr1.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0032_chr1.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0033_chr1.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0034_chr1.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0035_chr1.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0036_chr1.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0037_chr1.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0038_chr1.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0039_chr1.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0040_chr1.vcf.gz \
-d -R "0|0" \
> /lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/combined/20hrlr_rn6_chr1.vcf

bgzip -c 20hrlr_rn6_chr1.vcf > 20hrlr_rn6_chr1.vcf.gz

# do the weight cut for bilge's file
bcftools annotate -x INFO,^FORMAT/GT -r chr1 /nics/c/home/hanyou/UM/projects/rat_genome_03/vcf/Rats10XG.rn6.gVCFpool.6nt.chrs1-20.M.X.Y.sorted2.vcf.gz -o /lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/combined/GATK_20hrlr_rn6_chr1.vcf.gz -O z

# filter for SNPs
cd /lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/combined/
#bcftools view --types snps -o GATK_20hrlr_rn6_chr1_snps.vcf.gz -O z GATK_20hrlr_rn6_chr1.vcf.gz
bcftools view --exclude-types indels,mnps,other -o 20hrlr_rn6_chr1_snps.vcf.gz -O z 20hrlr_rn6_chr1.vcf.gz

# compare the two VCF files
vcf-compare 20hrlr_rn6_chr1.vcf.gz GATK_20hrlr_rn6_chr1.vcf.gz > compare.txt 
vcf-compare 20hrlr_rn6_chr1_snps.vcf.gz GATK_20hrlr_rn6_chr1_snps.vcf.gz > compare_no_indel.txt

# filter for minor allele frequency of 0.051
bcftools view -q 0.051:minor -O z 20hrlr_rn6_chr1_snps.vcf.gz -o 20hrlr_rn6_chr1_snps_maf.vcf.gz

# filter for phased only
bcftools view -q 0.051:minor --phased -O z 20hrlr_rn6_chr1_snps.vcf.gz -o 20hrlr_rn6_chr1_snps_maf_pha.vcf.gz

# filter for 10X filters
bcftools view -q 0.051:minor --phased -f .,PASS -O z 20hrlr_rn6_chr1_snps.vcf.gz -o 20hrlr_rn6_chr1_snps_maf_pha_f.vcf.gz

# Filter for QUAL 30
bcftools view -q 0.051:minor --phased -f .,PASS -i '%QUAL>30' -O z 20hrlr_rn6_chr1_snps.vcf.gz -o 20hrlr_rn6_chr1_snps_maf_pha_f_Q30.vcf.gz

# Filter for DP 10
bcftools view -q 0.051:minor --phased -f .,PASS -i 'FMT/DP>10' -O z 20hrlr_rn6_chr1_snps_maf_pha_f_Q30.vcf.gz -o 20hrlr_rn6_chr1_snps_maf_pha_f_Q30_DP10.vcf.gz

# make download file
bcftools view -H 20hrlr_rn6_chr1_snps_maf_pha_f.vcf.gz | grep -A 700 -P "\t8001461" > chr1_8mb+700.txt
bcftools view -h 20hrlr_rn6_chr1_snps_maf_pha_f.vcf.gz > chr1_8mb+700.vcf 
cat chr1_8mb+700.txt >> chr1_8mb+700.vcf

bcftools view -H -r chr1:8000000-108000000 20hrlr_rn6_chr1_snps_maf_pha_f_Q30_DP10_thin10k.recode.vcf.gz | head -700 > temp.txt
bcftools view -h 20hrlr_rn6_chr1_snps_maf_pha_f_Q30_DP10_thin10k.recode.vcf.gz > temp.vcf
cat temp.txt >> temp.vcf

# VCFtools, thin
vcftools --gzvcf 20hrlr_rn6_chr1_snps_maf_pha_f_Q30_DP10.vcf.gz --out 20hrlr_rn6_chr1_snps_maf_pha_f_Q30_DP10_thin10k --thin 10000 --recode
bgzip -c 20hrlr_rn6_chr1_snps_maf_pha_f_Q30_DP10_thin10k.recode.vcf > 20hrlr_rn6_chr1_snps_maf_pha_f_Q30_DP10_thin10k.recode.vcf.gz
tabix -p vcf 20hrlr_rn6_chr1_snps_maf_pha_f_Q30_DP10_thin10k.recode.vcf.gz


