#!/bin/bash
cd /lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/haploview


vcftools --vcf /lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/20hrlr_rn6_00.vcf --plink --chr chr1 --from-bp 8000000 --to-bp 8100000 --out t

vcftools --vcf temp01.vcf --plink --chr chr1 --from-bp 8000000 --to-bp 8100000 --out t

plink --file t --recode HV --snps-only just-acgt --out t


