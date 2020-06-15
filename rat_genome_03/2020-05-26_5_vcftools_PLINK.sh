#!/bin/bash
cd /lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/vcf/

vcftools --vcf /lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/20hrlr_rn6.vcf --plink --chr chr1 --from-bp 1 --to-bp 500000 --out t
