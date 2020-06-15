#!/bin/bash
/lustre/haven/proj/UTHSC0013/UM/apps/gatk_4.1.7.0/gatk-4.1.7.0/gatk CombineGVCFs \
-R /lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/references/refdata-rn6.fa.concatenate/fasta/genome.fa \
--variant /lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0021/5739-JL-0021.vcf \
--variant /lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0022/5739-JL-0022.vcf \
-O /lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/20hrlr_rn6.vcf



