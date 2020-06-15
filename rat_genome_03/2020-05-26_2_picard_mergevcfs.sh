#!/bin/bash
java -jar /lustre/haven/proj/UTHSC0013/UM/apps/picard_2.20.2/picard.jar MergeVcfs \
I=/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0021/5739-JL-0021.vcf \
I=/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0022/5739-JL-0022.vcf \
O=/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/20hrlr_rn6.vcf

