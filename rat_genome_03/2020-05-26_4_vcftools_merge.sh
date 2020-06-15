#!/bin/bash

#cd
#. .bashrc
#conda activate hanyou



# run after previous has been finished
for i in {21..40}
do
cd /lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/"5739-JL-00"$i/
bgzip -c "5739-JL-00"$i".vcf" > ../"5739-JL-00"$i".vcf.gz" &
done

# run after previous has been finished
for i in {21..40}
do
cd /lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs
tabix -p vcf "5739-JL-00"$i".vcf.gz" &
done

# run after previous has been finished
vcf-merge \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0021.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0022.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0023.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0024.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0025.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0026.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0027.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0028.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0029.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0030.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0031.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0032.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0033.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0034.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0035.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0036.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0037.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0038.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0039.vcf.gz \
/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/5739-JL-0040.vcf.gz \
> /lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/20hrlr_rn6.vcf

